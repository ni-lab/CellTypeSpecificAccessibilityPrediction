#!/usr/bin/env python

from optparse import OptionParser
import os
import json
from tqdm import tqdm
import pysam
import numpy as np
import pandas as pd

from basenji import seqnn
from basenji import dna_io

import tensorflow as tf
tf.compat.v1.flags.DEFINE_string('f','','')

'''
basenji_attributions.py
Compute gradient attributions for a set of input sequences.
'''

################################################################################
# main
################################################################################
def main():
  usage = 'usage: %prog [options] <model_dir> <bed_file>'
  parser = OptionParser(usage)
  parser.add_option('-f', dest='genome_fasta',
      default='%s/data/hg19.fa' % os.environ['BASENJIDIR'],
      help='Genome FASTA for sequences [Default: %default]')
  parser.add_option('-c', dest='cluster_name',default=None)
  parser.add_option('-n', dest='n', default=1000, type='int')
  parser.add_option('--random_state', dest='random_state',
                    default=42, type='int')
  parser.add_option('--rc', dest='rc',
      default=False, action='store_true',
      help='Average forward and reverse complement predictions [Default: %default]')
  parser.add_option('--shifts', dest='shifts',
      default='0', type='str',
      help='Ensemble prediction shifts [Default: %default]')
  (options, args) = parser.parse_args()

  if len(args) == 2:
    model_dir = args[0]
    bed_file = args[1]
  else:
    parser.error('Must provide parameters and model files and QTL VCF file')

  fasta_open = pysam.FastaFile(options.genome_fasta)

  # load model
  with open(f"{model_dir}/params.json") as params_open:
    params = json.load(params_open)
  params_model = params['model']
  options.shifts = [int(shift) for shift in options.shifts.split(',')]

  seqnn_model = seqnn.SeqNN(params_model)
  seqnn_model.restore(f"{model_dir}/model_best.h5")
  target_slice = np.arange(seqnn_model.model.output_shape[-1])
  seqnn_model.build_slice(target_slice)
  seqnn_model.build_ensemble(options.rc, options.shifts)
  model_ensemble = seqnn_model.ensemble

  cluster_seqs = pd.read_csv(bed_file, sep="\t", names=["chr", "start", "end", "peak_i"])
  if options.n > 0: 
    cluster_seqs = cluster_seqs.sample(n=options.n, replace=False, random_state=options.random_state) # sample 1,000 sequences from each peak cluster
  cluster_seqs.index = np.arange(len(cluster_seqs))
    
  # Get attributions using input * gradient
  attributions = np.zeros((len(cluster_seqs), 1344, 4, model_ensemble.output_shape[-1]), dtype=np.float32)
  one_hot_seqs = np.zeros((len(cluster_seqs), 1344, 4), dtype=int)

  for i, seq_row in tqdm(cluster_seqs.iterrows(), total=len(cluster_seqs)):
    seq = fasta_open.fetch(seq_row["chr"], seq_row["start"]-422, seq_row["end"]+422)
    seq_1hot = dna_io.dna_1hot(seq)
    
    attributions[i,:,:,:] = compute_gradient(model_ensemble, seq_1hot)
    one_hot_seqs[i,:,:] = seq_1hot
 
  inputxgrad = np.multiply(np.expand_dims(one_hot_seqs, axis=3), attributions).astype(np.float32)
  one_hot_seqs = one_hot_seqs.transpose(0, 2, 1)
  inputxgrad = inputxgrad.transpose(0, 2, 1, 3)
  np.save(f"{model_dir}/{options.cluster_name}_gradients.npy", attributions)
  np.save(f"{model_dir}/{options.cluster_name}_one_hot_seqs.npy", one_hot_seqs)
  # np.save(f"{model_dir}/{options.cluster_name}_inputxgradient.npy", inputxgrad)

  for ti in range(inputxgrad.shape[-1]):
    np.save(f"{model_dir}/{options.cluster_name}_inputxgradient_ti_{ti}.npy", inputxgrad[:,:,:,ti])
    
def compute_gradient(model, seq_1hot):
    # verify tensor shape
    seq_1hot = tf.convert_to_tensor(seq_1hot, dtype=tf.float32)
    if len(seq_1hot.shape) < 3:
        seq_1hot = tf.expand_dims(seq_1hot, axis=0)

    # sequence input
    sequence = tf.keras.Input(shape=(seq_1hot.shape[1], 4), name='sequence')

    # predict
    predictions = model(sequence)
    
    # slice
    target_slice = np.arange(model.output_shape[-1])
    predictions_slice = tf.gather(predictions, target_slice, axis=-1)
    predictions_slice = tf.gather(predictions_slice, 0, axis=-2)
    
    # replace model
    model_batch = tf.keras.Model(inputs=sequence, outputs=predictions_slice)
    
    with tf.GradientTape(persistent=True) as tape:
        tape.reset()
        tape.watch(seq_1hot)
        # predict
        preds = model_batch(seq_1hot, training=False)

        # sum across positions
        preds = tf.reduce_sum(preds, axis=-2)

    # compute jacobian
    grads = tape.jacobian(preds, seq_1hot, experimental_use_pfor=False)
    grads = tf.squeeze(grads, [1])
    grads = tf.transpose(grads, [1,2,0])
    
    return grads.numpy()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
  main()
