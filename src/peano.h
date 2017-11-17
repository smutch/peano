long long peano_hilbert_key(int x, int y, int z, int bits);
void peano_hilbert_keys(int* x, int* y, int* z, int n_points, int bits, long long* keys);

void peano_hilbert_key_inverse(long long key, int bits, int *x, int *y, int *z);
void peano_hilbert_keys_inverse(long long* keys, int n_keys, int bits, int *x, int *y, int *z);
