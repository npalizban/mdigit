/*************************
 *
 * This file implements fixed point and above 64-bit arithmatic functions.
 *
 * Note 1: Granularity (each digit) is 32 bit:
 *
 * Note 2: Implementation is based on big-endian
 *
 *************************/

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include "mdigit.h"

uint64_t __udivmoddi4(uint64_t u, uint64_t v, uint64_t* r);

/*
 * in place two's complement of m digit number
 *
 * param in/out     a       input number
 * param in         m       number of digits of a
 *
 */
void negate_m(uint32_t* a, uint32_t m)
{
    for (uint32_t i = 0; i < m; i++)
    {
        a[i] = ~a[i];
    }

    uint32_t k = m - 1;
    a[k]++;

    while ((a[k] == 0) && k > 0)
    {
        k--;
        a[k]++;
    }
}


/*
 * in place left shift m digit number.
 * note: left shift is defined as shift from LSB to MSB.
 *
 * param in/out     a           input number
 * param in         m           number of digits of a
 * param in         shift       number of bits to shift
 *
 */
void lshift_m(uint32_t* a, uint32_t m, uint32_t shift)
{
    if (shift > (m << 5))
    {
        shift = m << 5;
    }

    int32_t k = shift >> 5; // number of full digits to shift
    shift &= 0x1F; // shift % 32

    for (int32_t i = 0; i < (int32_t) m - k - 1; i++)
    {
        a[i] = (a[i + k] << shift) | (((uint64_t) a[i + k + 1]) >> (32 - shift));
    }

    if (m > k)
    {
        a[m - k - 1] = a[m - 1] << shift;
    }

    for (int32_t i = m - k; i < m ; i++)
    {
        a[i] = 0;
    }
}


/*
 * in place right shift of m digit number
 *
 * param in/out     a
 * param in         m       number of digits in a
 * param in         shift   number of bits to shift
 * param in         sext    sign extend flag. If false append zero, if true sign extend.
 */
void rshift_m(uint32_t* a, uint32_t m, uint32_t shift, uint32_t sext)
{
    if (shift > m << 5)
    {
        shift = m << 5;
    }
    int32_t k = shift >> 5; // number of full digits to shift
    shift &= 0x1F; // shift % 32

    for (uint32_t i = m - 1; i > k; i--)
    {
        a[i] = (((uint64_t) a[i - k]) >> shift) | ((uint64_t) (a[i - k - 1]) << (32 - shift));
    }

    if (sext && (a[0] >> 31))
    {
        if (k < m)
        {
            a[k] = ((int32_t) a[0]) >> shift;
        }

        for (int32_t i = 0; i < k; i++)
        {
            a[i] = 0xffffffff;
        }
    }
    else
    {
        if (k < m)
        {
            a[k] = a[0] >> shift;
        }

        for (int32_t i = 0; i < k; i++)
        {
            a[i] = 0x00;
        }
    }
}

/*
 * add m digit numbers and store it in m digit result.
 *
 * param in     a       operand 1
 * param in     b       operand 2
 * param in     m       number of digits of a and b
 * param out    res     result of a + b
 *
 */
void add_m(uint32_t* a, uint32_t* b, uint32_t m, uint32_t* res)
{
    uint32_t carry = 0;
    for (int32_t i = m - 1; i >= 0; i--)
    {
        uint64_t add = (((uint64_t) a[i]) + b[i]) + carry;

        res[i] = add;
        carry = add >> 32;
    }
}

/*
 * add m digit numbers and store it in m digit result.
 *
 * param in     a       operand 1
 * param in     b       operand 2
 * param in     c       operand 3
 * param in     m       number of digits of a and b and c
 * param out    res     result of a + b + c
 *
 */
void add3_m(uint32_t* a, uint32_t* b, uint32_t* c, uint32_t m, uint32_t* res)
{
    uint32_t carry = 0;
    for (int32_t i = m - 1; i >= 0; i--)
    {
        uint64_t add = ((((uint64_t) a[i]) + b[i]) + c[i]) + carry;
        res[i] = add;
        carry = add >> 32;
    }
}


/*
 * subtract signed m digit numbers and store in m digit result
 *
 * param in     a       operand 1
 * param in     b       operand 2
 * param in     m       number of digits of a and b
 * param out    res     result of a - b
 */
void sub_m(uint32_t* a, uint32_t* b, uint32_t m, uint32_t* res)
{
    uint32_t borrow = 0;
    for (int32_t i = m - 1; i >= 0; i--)
    {
        uint64_t sub = (((uint64_t) a[i]) - b[i]) - borrow;
        res[i] = sub;
        borrow = (sub >> 63);
    }
}

/*
 * unsigned multiply m digit by n digit and
 * store it in m + n digit result
 *
 * param in     u       operand 1
 * param in     v       operand 2
 * param in     m       number of digits of u
 * param in     n       number of digits of v
 * param in     sign    if 0, usigned multiply, if 1 only u is signed, if 2 u and v are signed,
 * param out    res     m + n digit u * v
 *
 * if an operand is signed and negative, we have:
 * u * v = -(-u * v) = -((2^m - u) * v = u * v - 2^m * v
 */
void mul_m_n(uint32_t* u, uint32_t* v, uint32_t m, uint32_t n, uint32_t sign, uint32_t* res)
{
    uint64_t sum;
    uint32_t carry = 0;
    for (int32_t i = n - 1; i >= 0; i--)
    {
        sum = ((uint64_t) u[m - 1]) * v[i] + carry;
        res[i + m] = sum;
        carry = sum >> 32;
    }
    res[m - 1] = carry;

    for (int32_t i = m - 2; i >= 0; i--)
    {
        carry = 0;
        for (int32_t j = n - 1; j >= 0; j--)
        {
            sum = (((uint64_t) u[i]) * v[j] + carry) + res[i + j + 1];
            res[i + j + 1] = sum;
            carry = sum >> 32;
        }
        res[i] = carry;
    }

    if (sign)
    {
        if ((int32_t) u[0] < 0)
        {
            uint64_t sub;
            uint32_t b = 0;
            for (int32_t i = n - 1; i >= 0; i--)
            {
                sub = (((uint64_t) res[i]) - v[i]) - b;
                res[i] = sub;
                b = sub >> 63;
            }
        }
        if (sign > 1 && ((int32_t) v[0]) < 0)
        {
            uint64_t sub;
            uint32_t b = 0;
            for (int32_t i = m - 1; i >= 0; i--)
            {
                sub = (((uint64_t) res[i]) - u[i]) - b;
                res[i] = sub;
                b = sub >> 63;
            }
        }
    }
}

/*
 * multiply m by n digits and add the result to m + n digit result.
 *
 * param in         u       operand 1
 * param in         v       operand 2
 * param in         m       number of digits of u
 * param in         n       number of digits of v
 * param in         sign    if 0, usigned multiply, if 1 only u is signed, if 2 u and v are signed,
 * param in/out     res     m + n result. res = res + (u * v)
 *
 */
void mad_m_n(uint32_t* u, uint32_t* v, uint32_t m, uint32_t n, uint32_t sign, uint32_t* res)
{
    for (int32_t i = n - 1; i >= 0; i--)
    {
        uint32_t carry = 0;
        for (int32_t j = m - 1; j >= 0; j--)
        {
            uint64_t add = ((uint64_t) u[j]) * v[i] +
                           ((uint64_t) carry) +
                           ((uint64_t) res[i + j + 1]);

            res[i + j + 1] = add;
            carry = add >> 32;
        }
        res[i] += carry;
    }

    if (sign)
    {
        if (((int32_t) u[0]) < 0)
        {
            uint64_t sub;
            uint32_t b = 0;
            for (int32_t i = n - 1; i >= 0; i--)
            {
                sub = (((uint64_t) res[i]) - v[i]) - b;
                res[i] = sub;
                b = sub >> 63;
            }
        }
        if (sign > 1 && ((int32_t) v[0]) < 0)
        {
            uint64_t sub;
            uint32_t b = 0;
            for (int32_t i = m - 1; i >= 0; i--)
            {
                sub = (((uint64_t) res[i]) - u[i]) - b;
                res[i] = sub;
                b = sub >> 63;
            }
        }
    }
}

/*
 *  Divid unsigned m digit dividend by n digit divisor into
 *  m digit quotient and n digit reminder in base 2^32 with
 *  Knuth algorithm.
 *
 *  Implementation is based on Hacker's delight book.
 *
 *  notes:
 *
 *      1. n, m <= DIV_KNUTH_MAX_DIGITS
 *      2. if divisor is zero, function returns maximum value for
 *         reminder and quotient.
 *
 * param in     u           dividend
 * param in     v           divisor
 * param in     m           number of digits of u
 * param in     n           number of digits of v
 * param out    q           quotient
 * param out    r           reminder
 *
 */
void div_m_n(uint32_t* u,  uint32_t* v, uint32_t m, uint32_t n, uint32_t* q, uint32_t* r)
{
    // store initial value of n then remove
    // most significant zero digits from divisor.
    uint32_t n_copy = n;
    while (v[0] == 0 && n > 0)
    {
        n--;
        v++;
    }
    if (n == 1)
    {
        /* if divisor is a single digit, don't use Knuth, instead, expand the dividend. */
        uint32_t rh = 0;
        for (uint32_t i = 0; i < m; i++)
        {
            uint64_t un = (((uint64_t) rh) << 32) + u[i];
            q[i] = un / v[0];
            rh = un - ((uint64_t) q[i]) * v[0];
        }

        if (r != NULL)
        {
            for (uint32_t i = 0; i < n_copy - n; i++)
            {
                *r = 0;
                r++;
            }

            *r = rh;
        }
        return;
    }
    else if (n > m)
    {
        for (uint32_t i = 0; i < m; i++)
        {
            q[i] = 0x00UL;
        }
        if (r != NULL)
        {
            for (uint32_t i = 0; i < n_copy - m; i++)
            {
                r[i] = 0x00UL;
            }
            for (uint32_t i = 0; i < m; i++)
            {
                r[n_copy - m + i] = u[i];
            }
        }
        return;
    }
    else if (n == 0)
    {
        /* division by zero, return maximum */
        for (uint32_t i = 0; i < m; i++)
        {
            q[i] = 0xFFFFFFFFUL;
        }
        if (r != NULL)
        {
            for (uint32_t i = 0; i < n_copy; i++)
            {
                r[i] = 0xFFFFFFFFUL;
            }
        }
        return;
    }


    /*
     * normalize
     */
    uint32_t vn[DIV_KNUTH_MAX_DIGITS];
    uint32_t un[DIV_KNUTH_MAX_DIGITS + 1]; // extra one digit for normalization overflow

    // number of leading zeros.
    // max nlz is 31 since v[0] > 0;
    uint32_t nlz = __builtin_clz(v[0]);

    // number of leading zeros complement.
    // note: nlz_cmp can be 32 if nlz is zero therefore
    // when shifting by nlz_cmp need to cast to 64 bit first
    // max nlz_cmp = 32
    uint32_t nlz_cmp = 32 - nlz;


    // vn = v << nlz
    for (uint32_t i = 0; i < n - 1; i++)
    {
        vn[i] = (v[i] << nlz) | (((uint64_t )v[i + 1]) >> nlz_cmp);
    }
    vn[n - 1] = v[n - 1] << nlz;

    // un = u << nlz
    // un is m + 1 digits because of potential overflow
    un[0] = ((uint64_t) u[0]) >> nlz_cmp;
    for (uint32_t i = 0; i < m - 1; i++)
    {
        un[i + 1] = (u[i] << nlz) | (((uint64_t)u[i + 1]) >> nlz_cmp);
    }
    un[m] = u[m - 1] << nlz;

    // quotient has m - n + 1 digits
    // therefore, the most significant m - (m - n + 1) = n - 1 digits are zero
    for (uint32_t i = 0; i < n - 1; i++)
    {
        q[i] = 0;
    }

    /*
     * long division. each iteration
     * divides n + 1 by n digit
     */
    for (uint32_t j = 0; j <= m - n; j++)
    {

        // un_p is n + 1 digit portion of dividend.
        // In each iteration, the n most significant digits of un_p
        // is the remainder of previous iteration. Therefore, un_p[0..n-1] < vn[0..n-1].
        // un_p[n] is the new digit peeled from un.
        uint32_t* un_p = &un[j];

        // Estimage the quotient
        uint64_t rh;
        uint64_t qh = __udivmoddi4(((uint64_t) un_p[0] << 32) + un_p[1], vn[0], &rh);

        const uint64_t b = (1ULL << 32);

        // Below loop corrects the estimate by looking at divisor second digit.
        // After this loop, qh is guranteed to be exact or at most off by 1
        // and it is guranteed to be smaller than b.
        // max iterations of this loop is 2
        while ((qh >= b) || ((qh * vn[1]) > (b * rh + un_p[2])))
        {
            qh = qh - 1;
            rh = rh + vn[0];
            if (rh >= b)
            {
                break;
            }
        }

        // Multiply and subtract.
        // un_p -= qh * vn
        int64_t k = 0;
        int64_t t;
        for (uint32_t i = n; i > 0; i--)
        {
            uint64_t p = (uint64_t) qh * vn[i - 1];

            t = ((uint64_t) un_p[i]) - k - (p & 0xFFFFFFFF);
            un_p[i] = t;
            k = (p >> 32) - (t >> 32);
        }

        t = un_p[0] - k;
        un_p[0] = t;

        // If we subtracted too much, add back
        // un_p += vn
        if (t < 0)
        {
            qh = qh - 1;
            k = 0;
            for (uint32_t i = n; i > 0; i--)
            {
                t = ((uint64_t) un_p[i]) + vn[i - 1] + k;
                un_p[i] = t;
                k = t >> 32;
            }
            un_p[0] = un_p[0] + k;
        }

        q[j + n - 1] = qh;
    }

    /*
     * unnormalize the reminder if wanted
     */
    if (r != NULL)
    {
        for (uint32_t i = 0; i < n_copy - n; i++)
        {
            r[i] = 0;
        }

        for (uint32_t i = 0; i < n - 1; i++)
        {
            r[n_copy - i - 1] = (((uint64_t) un[m - i]) >> nlz) |
                                ((uint64_t) (un[m - i - 1]) << nlz_cmp);
        }

        r[n_copy - n] = ((uint64_t) un[m - n + 1]) >> nlz;
    }
}

/*
 * divide m digit by n digit and store the result
 * in k integer part + f fractional part digit quotient
 *
 * requires 1) m + f <= DIV_KNUTH_MAX_DIGITS
 *          2) k <= m
 *
 * param in     u       dividend
 * param in     v       divisor
 * param in     m       number of digits of u
 * param in     n       number of digits of v
 * param in     k       number of integer digits in q
 * param in     f       number of fractional digits in q
 * param out    q       quotient with k real digits and f fractional digits
 */
void div_m_n_fract_k_f(uint32_t* u, uint32_t* v, uint32_t m, uint32_t n, uint32_t k, uint32_t f, uint32_t* q)
{
    // s = u << (32 * f)
    uint32_t s[DIV_KNUTH_MAX_DIGITS];
    for (uint32_t i = 0; i < m; i++)
    {
        s[i] = u[i];
    }
    for (uint32_t i = m; i < m + f; i++)
    {
        s[i] = 0;
    }

    div_m_n(s, v, m + f, n, s, NULL);

    for (int32_t i = 0; i < k + f; i++)
    {
        q[k + f - i - 1] = s[m + f - i - 1];
    }
}
