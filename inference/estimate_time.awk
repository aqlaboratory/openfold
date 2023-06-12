#!/bin/awk -f
#
# Copyright 2023 RIKEN & Fujitsu Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

BEGIN{
    coeff_a = 0.0023
    coeff_b = 0.4875
    coeff_c = 35.636
    max_seq_len = 2000

    n_redundant_seq = 0
}

/^[^>].*/{
    seqs[$0] = 0
    n_redundant_seq += 1
}
END{
    n_seq = 0
    n_acc_seq = 0
    n_unk_seq = 0
    acc_est_time = 0

    for (seq in seqs) {
	seq_len = length(seq)

	n_unk = split(seq, tmp, "X") - 1
    
	if (seq_len <= max_seq_len && n_unk == 0) {
	    est_time = coeff_a * seq_len * seq_len + coeff_b * seq_len + coeff_c
	    acc_est_time += est_time
	    n_acc_seq += 1
	}
	if (n_unk > 0) {
	    n_unk_seq += 1
	}
	n_seq += 1
    }

    printf("Time estimation function: %f x^2 + %f x + %f\n", coeff_a, coeff_b, coeff_c)
    printf("Number of seq: %d\n", n_redundant_seq)
    printf("Number of unique seq: %d\n", n_seq)
    printf("Number of unknown seq: %d\n", n_unk_seq)
    printf("Number of time accumulated seq: %d\n", n_acc_seq)
    printf("Estimated time [hours]: %f\n", acc_est_time/3600)
}
