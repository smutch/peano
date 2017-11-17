#include <criterion/criterion.h>
#include <criterion/logging.h>
#include <criterion/parameterized.h>
#include "../src/peano.h"

Test(single, zero)
{
    long long key = peano_hilbert_key(0, 0, 0, 8);
    cr_assert_eq(key, 0ll);

    int x, y, z;
    peano_hilbert_key_inverse(key, 8, &x, &y, &z);
    cr_assert_eq(x, 0);
    cr_assert_eq(y, 0);
    cr_assert_eq(z, 0);
}

Test(single, forward_back)
{
    int x = 7;
    int y = 126;
    int z = 255;
    long long key = peano_hilbert_key(x, y, z, 8);
    cr_assert_eq(key, 14680266ll);

    int xr, yr, zr;
    peano_hilbert_key_inverse(key, 8, &xr, &yr, &zr);
    cr_assert_eq(xr, x);
    cr_assert_eq(yr, y);
    cr_assert_eq(zr, z);
}

Test(single, fail_low_bits, .signal = SIGABRT)
{
    peano_hilbert_key(0, 0, 28, 3);
}

Test(single, max_bits)
{
    peano_hilbert_key(0, 0, 26, 3);
}

Test(array, dummy)
{
    const unsigned n_points = 5;
    long long keys[n_points];
    int x[n_points], y[n_points], z[n_points];

    for (unsigned ii = 0; ii < n_points; ii++) {
        x[ii] = 0;
        y[ii] = 0;
        z[ii] = 0;
    }

    peano_hilbert_keys(x, y, z, (int)n_points, 8, keys);

    for (unsigned ii = 0; ii < n_points; ii++) {
        cr_log_info("x[%d]=%d, y[%d]=%d, z[%d]=%d => key[%d]=%lld\n", ii, x[ii], ii, y[ii], ii, z[ii], ii, keys[ii]);
    }
}

ParameterizedTestParameters(array, zero_forward)
{
    const unsigned n_points = 5;
    static long long keys[n_points];  // N.B. use of static here so as available outside
    int x[n_points], y[n_points], z[n_points];

    for (unsigned ii = 0; ii < n_points; ii++) {
        x[ii] = 0;
        y[ii] = 0;
        z[ii] = 0;
    }

    peano_hilbert_keys(x, y, z, (int)n_points, 8, keys);

    return cr_make_param_array(long long, keys, (int)n_points);
}

ParameterizedTest(long long* key, array, zero_forward)
{
    cr_log_info("key = %lld\n", *key);
    cr_assert_eq(*key, 0ll);
}
