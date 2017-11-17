#include <assert.h>
#include <math.h>

/*
 * This code was taken and minimally adapted from the GADGET-2 hydrodynamic
 * cosmological simulations code.
 * (https://wwwmpa.mpa-garching.mpg.de/gadget/html/index.html)
 *
 * Author: Volker Springel
 */

static const int quadrants[24][2][2][2] = {
    /* rotx=0, roty=0-3 */
    {{{0, 7}, {1, 6}}, {{3, 4}, {2, 5}}},
    {{{7, 4}, {6, 5}}, {{0, 3}, {1, 2}}},
    {{{4, 3}, {5, 2}}, {{7, 0}, {6, 1}}},
    {{{3, 0}, {2, 1}}, {{4, 7}, {5, 6}}},
    /* rotx=1, roty=0-3 */
    {{{1, 0}, {6, 7}}, {{2, 3}, {5, 4}}},
    {{{0, 3}, {7, 4}}, {{1, 2}, {6, 5}}},
    {{{3, 2}, {4, 5}}, {{0, 1}, {7, 6}}},
    {{{2, 1}, {5, 6}}, {{3, 0}, {4, 7}}},
    /* rotx=2, roty=0-3 */
    {{{6, 1}, {7, 0}}, {{5, 2}, {4, 3}}},
    {{{1, 2}, {0, 3}}, {{6, 5}, {7, 4}}},
    {{{2, 5}, {3, 4}}, {{1, 6}, {0, 7}}},
    {{{5, 6}, {4, 7}}, {{2, 1}, {3, 0}}},
    /* rotx=3, roty=0-3 */
    {{{7, 6}, {0, 1}}, {{4, 5}, {3, 2}}},
    {{{6, 5}, {1, 2}}, {{7, 4}, {0, 3}}},
    {{{5, 4}, {2, 3}}, {{6, 7}, {1, 0}}},
    {{{4, 7}, {3, 0}}, {{5, 6}, {2, 1}}},
    /* rotx=4, roty=0-3 */
    {{{6, 7}, {5, 4}}, {{1, 0}, {2, 3}}},
    {{{7, 0}, {4, 3}}, {{6, 1}, {5, 2}}},
    {{{0, 1}, {3, 2}}, {{7, 6}, {4, 5}}},
    {{{1, 6}, {2, 5}}, {{0, 7}, {3, 4}}},
    /* rotx=5, roty=0-3 */
    {{{2, 3}, {1, 0}}, {{5, 4}, {6, 7}}},
    {{{3, 4}, {0, 7}}, {{2, 5}, {1, 6}}},
    {{{4, 5}, {7, 6}}, {{3, 2}, {0, 1}}},
    {{{5, 2}, {6, 1}}, {{4, 3}, {7, 0}}}
};


static const int rotxmap_table[24] = { 4, 5, 6, 7, 8, 9, 10, 11,
    12, 13, 14, 15, 0, 1, 2, 3, 17, 18, 19, 16, 23, 20, 21, 22
};

static const int rotymap_table[24] = { 1, 2, 3, 0, 16, 17, 18, 19,
    11, 8, 9, 10, 22, 23, 20, 21, 14, 15, 12, 13, 4, 5, 6, 7
};

static const int rotx_table[8] = { 3, 0, 0, 2, 2, 0, 0, 1 };
static const int roty_table[8] = { 0, 1, 1, 2, 2, 3, 3, 0 };

static const int sense_table[8] = { -1, -1, -1, +1, +1, -1, -1, -1 };

static int flag_quadrants_inverse = 1;
static char quadrants_inverse_x[24][8];
static char quadrants_inverse_y[24][8];
static char quadrants_inverse_z[24][8];


long long peano_hilbert_key(int x, int y, int z, int bits)
{
    int max = pow(3, bits);
    assert(x < max);
    assert(y < max);
    assert(z < max);

    int mask = 1 << (bits - 1);
    long long key = 0;
    int rotation = 0;
    int sense = 1;

    for(int ii = 0; ii < bits; ii++, mask >>= 1)
    {
        const int bitx = (x & mask) ? 1 : 0;
        const int bity = (y & mask) ? 1 : 0;
        const int bitz = (z & mask) ? 1 : 0;

        const int quad = quadrants[rotation][bitx][bity][bitz];

        key <<= 3;
        key += (sense == 1) ? (quad) : (7 - quad);

        int rotx = rotx_table[quad];
        int roty = roty_table[quad];
        sense *= sense_table[quad];

        while(rotx > 0)
        {
            rotation = rotxmap_table[rotation];
            rotx--;
        }

        while(roty > 0)
        {
            rotation = rotymap_table[rotation];
            roty--;
        }
    }

    return key;
}


void peano_hilbert_key_inverse(long long key, int bits, int *x, int *y, int *z)
{

    if(flag_quadrants_inverse)
    {
        flag_quadrants_inverse = 0;
        for(int rotation = 0; rotation < 24; rotation++)
            for(int bitx = 0; bitx < 2; bitx++)
                for(int bity = 0; bity < 2; bity++)
                    for(int bitz = 0; bitz < 2; bitz++) {
                        int quad = quadrants[rotation][bitx][bity][bitz];
                        quadrants_inverse_x[rotation][quad] = bitx;
                        quadrants_inverse_y[rotation][quad] = bity;
                        quadrants_inverse_z[rotation][quad] = bitz;
                    }
    }

    int shift = 3 * (bits - 1);
    int mask = 7 << shift;

    int rotation = 0;
    char sense = 1;

    *x = *y = *z = 0;

    for(int ii = 0; ii < bits; ii++, mask >>= 3, shift -= 3)
    {
        int keypart = (key & mask) >> shift;

        int quad = (sense == 1) ? (keypart) : (7 - keypart);

        *x = (*x << 1) + quadrants_inverse_x[rotation][quad];
        *y = (*y << 1) + quadrants_inverse_y[rotation][quad];
        *z = (*z << 1) + quadrants_inverse_z[rotation][quad];

        int rotx = rotx_table[quad];
        int roty = roty_table[quad];
        sense *= sense_table[quad];

        while(rotx > 0)
        {
            rotation = rotxmap_table[rotation];
            rotx--;
        }

        while(roty > 0)
        {
            rotation = rotymap_table[rotation];
            roty--;
        }
    }
}
