#include <criterion/criterion.h>
#include <criterion/logging.h>
#include "../src/peano.h"

Test(basic, zero)
{
    long long key = peano_hilbert_key(0, 0, 0, 8);
    cr_assert_eq(key, 0ll);

    int x, y, z;
    peano_hilbert_key_inverse(key, 8, &x, &y, &z);
    cr_assert_eq(x, 0);
    cr_assert_eq(y, 0);
    cr_assert_eq(z, 0);
}

Test(basic, forward_back)
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

Test(basic, fail_low_bits, .signal = SIGABRT)
{
    peano_hilbert_key(0, 0, 28, 3);
}

Test(basic, max_bits)
{
    peano_hilbert_key(0, 0, 26, 3);
}
