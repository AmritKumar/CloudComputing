/*
 *  test.h
 *  integer-fhe
 *
 *  Created by Henning Perl on 17.12.10.
 *
 */

#pragma once
#ifndef TEST_H
#define TEST_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "integer-fhe.h"

#define RUNS 60
#define KEYRUNS 10

void test_encryt_decrypt();

void test_halfadd();

void test_fulladd();

void test_recrypt();

void test_homomorphic();

void test_fully_homomorphic();

void test_xor_bits();

void test_sum_bits();

void test_bit_majoritaire();

void test_suite();

#endif