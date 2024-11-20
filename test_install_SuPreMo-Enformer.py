#!/usr/bin/env python
# coding: utf-8
import random
from typing import Any
import numpy as np
import tensorflow.compat.v2 as tf
import tensorflow_hub as hub

def one_hot_encode(sequence: str,
                   alphabet: str = 'ACGT',
                   neutral_alphabet: str = 'N',
                   neutral_value: Any = 0,
                   dtype=np.float32) -> np.ndarray:
  """One-hot encode sequence."""
  def to_uint8(string):
    return np.frombuffer(string.encode('ascii'), dtype=np.uint8)
  hash_table = np.zeros((np.iinfo(np.uint8).max, len(alphabet)), dtype=dtype)
  hash_table[to_uint8(alphabet)] = np.eye(len(alphabet), dtype=dtype)
  hash_table[to_uint8(neutral_alphabet)] = neutral_value
  hash_table = hash_table.astype(dtype)
  return hash_table[to_uint8(sequence)]

#enformer_model = hub.load("https://kaggle.com/models/deepmind/enformer/frameworks/TensorFlow2/variations/enformer/versions/1").model
enformer_model = hub.load("/pollard/home/xzhang/.cache/kagglehub/models/deepmind/enformer/tensorFlow2/enformer/1").model
SEQ_LENGTH = 393_216
# SEQ_LENGTH = 114_688

# Input array [batch_size, SEQ_LENGTH, 4] one hot encoded in order 'ACGT'. The
# `one_hot_encode`. With N values being all zeros.
inputs = tf.zeros((1, SEQ_LENGTH, 4), dtype=tf.float32)

predictions = enformer_model.predict_on_batch(inputs)
predictions['human'].shape  # [batch_size, 896, 5313]
predictions['mouse'].shape  # [batch_size, 896, 1643]

sequences = []
for i in range(2):
    seq = ''.join([random.choice('ACGT') for _ in range(SEQ_LENGTH)])
    sequences.append(seq)

inputs = [np.expand_dims(one_hot_encode(s), 0).astype(np.float32) for s in sequences]
inputs = np.vstack(inputs)
predictions = enformer_model.predict_on_batch(inputs)
predictions['human'].shape  # [batch_size, 896, 5313]
predictions['mouse'].shape  # [batch_size, 896, 1643]

print(predictions['human'])
