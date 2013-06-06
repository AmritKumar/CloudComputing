/*
 *  test.c
 *  integer-fhe
 *
 *  Created by Henning Perl on 17.12.10.
 *
 */

#include "test.h"

#define ASSERT_HALFADD(__a,__b,__sum,__carry)		\
	fhe_halfadd(sum, carry, __a, __b, pk);			\
	assert(fhe_decrypt(sum, sk) == __sum);			\
	assert(fhe_decrypt(carry, sk) == __carry);

#define ASSERT_FULLADD(__a,__b,__c,__sum,__carry)	\
	fhe_fulladd(sum, carry, __a, __b, __c, pk);		\
	assert(fhe_decrypt(sum, sk) == __sum);			\
	assert(fhe_decrypt(carry, sk) == __carry);

#define ASSERT_HOMMUL(__a, __b, __check)			\
	fhe_mul(temp, __a, __b, pk);					\
	assert(fhe_decrypt(temp, sk) == __check);

#define ASSERT_HOMADD(__a, __b, __check)			\
	fhe_add(temp, __a, __b, pk);					\
	assert(fhe_decrypt(temp, sk) == __check);
void debug_test_bit_majoritaire();
void
test_suite()
{
	//test_fully_homomorphic();
	//test_homomorphic();
	//test_recrypt();
	//test_encryt_decrypt();
	//test_halfadd();
	//test_fulladd();
	test_xor_bits();
	//test_sum_bits();
	//test_bit_majoritaire();
	//debug_test_bit_majoritaire();
}

void debug_or(mpz_t res, mpz_t a, mpz_t b){
	mpz_t aux1, aux2;
	mpz_init(aux1);
	mpz_init(aux2);
	mpz_add(aux1, a, b);
	mpz_mul(aux2, a, b);
	mpz_add(res, aux1, aux2);
	mpz_clear(aux1);
	mpz_clear(aux2);
}

void or(mpz_t res, mpz_t a, mpz_t b, fhe_pk_t pk){
	mpz_t aux1, aux2;
	mpz_init(aux1);
	mpz_init(aux2);
	fhe_add(aux1, a, b, pk);
	fhe_mul(aux2, a, b, pk);
	fhe_add(res, aux1, aux2, pk);
	mpz_clear(aux1);
	mpz_clear(aux2);
}

void debug_not(mpz_t res, mpz_t a){
	mpz_t c1;
	mpz_init(c1);
	mpz_set_ui(c1, 1);
	mpz_sub(res, c1, a);
	mpz_clear(c1);
}
void not(mpz_t res, mpz_t a, fhe_pk_t pk){
	mpz_t aux, c1;
	mpz_init(aux);
	mpz_init(c1);
	fhe_encrypt(c1, pk, 1);
	fhe_add(res, a, c1, pk);
	mpz_clear(aux);
	mpz_clear(c1);
}
void test_xor_bits(){
	int  aux, i;
	mpz_t c0, c1;
	fmpz_poly_t poly_c;
	fmpz_poly_init(poly_c);
	mpz_init(c0);
	mpz_init(c1);
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	fhe_keygen(pk, sk);
	//m= 13;	
	for(int m=13; m < 27; m++){
		aux=m;
		i=0;
		printf("--------->%i\n", m % 2);
		fhe_encrypt(c0, pk, aux % 2);
		fmpz_poly_set_coeff_mpz ( poly_c , i , c0 );
		aux = aux >> 1;
		do {
			printf("--------->%i\n", aux % 2);
			fhe_encrypt(c0, pk, aux % 2);
			i=i+1;
			fmpz_poly_set_coeff_mpz ( poly_c , i , c0 );
			aux = aux >> 1;
		}while(aux != 0);
		fmpz_poly_get_coeff_mpz ( c0, poly_c, i);
		i=i-1;
		while(i>=0){
			fmpz_poly_get_coeff_mpz ( c1, poly_c, i);
			fhe_add(c0, c0, c1, pk);
			i=i-1;
		}
		aux = fhe_decrypt(c0, sk);
		printf("///////////////-----------!!!!!!  xor bits de %i est %i \n", m, aux);
	}
	mpz_clear(c0);
	mpz_clear(c1);	
	fmpz_poly_clear(poly_c);
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);
}
void debug_test_bit_majoritaire(){
	int a, res, aux, i;
	a=11;
	mpz_t c0;
	fmpz_poly_t poly_c;
	for(int s=19; s< 20; s++){
		a=s;
		i=0;
		mpz_init(c0);
		fmpz_poly_init( poly_c );
		printf("--------->%i\n", a % 2);
		mpz_set_ui(c0,  a % 2);
		fmpz_poly_set_coeff_mpz( poly_c , i , c0 );
		aux = a;	
		aux = aux >> 1;
		do {
			printf("--------->%i\n", aux % 2);
			mpz_set_ui(c0,  aux % 2);
			i=i+1;
			fmpz_poly_set_coeff_mpz ( poly_c , i , c0 );
			aux = aux >> 1;

		}while(aux != 0);
		mpz_t tmp1, tmp2, tmp3;
		mpz_init(tmp1);
		mpz_init(tmp2);
		mpz_init(tmp3);
		int k;
//première itération
		fmpz_poly_get_coeff_mpz ( tmp3 , poly_c , 0 );
		gmp_printf ("P(0)= %Zd\n", tmp3);
//
		for(int j=1; j<= i; j++){
			fmpz_poly_get_coeff_mpz ( tmp2 , poly_c , j );
			gmp_printf ("P( %i ) = %Zd\n", j, tmp2);
			mpz_mul(tmp3, tmp3, tmp2);
			gmp_printf (" %Zd\n", tmp3);
		}
//
		for(int j= i; j>=0; j--){
			fmpz_poly_get_coeff_mpz ( tmp1 , poly_c , j );
			gmp_printf ("P( %i ) = %Zd\n", j, tmp1);
			debug_not(tmp1, tmp1);
			gmp_printf ("notP( %i ) = %Zd\n", j, tmp1);

			k=j-1;
			while(k >=0){
				fmpz_poly_get_coeff_mpz ( tmp2 , poly_c , k );
				gmp_printf ("P( %i ) = %Zd\n", k, tmp2);
				mpz_mul(tmp1, tmp1, tmp2);
				gmp_printf ("!! %Zd\n", tmp1);
				k=k-1;
			}
			k= j+1;
			while(k <= i){
				fmpz_poly_get_coeff_mpz ( tmp2 , poly_c , k );
				gmp_printf ("P( %i ) = %Zd\n",k, tmp2);
				mpz_mul(tmp1, tmp1, tmp2);
				gmp_printf (" %Zd\n", tmp1);
				k=k+1;
			}
			debug_or(tmp3, tmp1, tmp3);
			gmp_printf (" %Zd\n", tmp3);
		}

		res = mpz_get_ui(tmp3);
		printf("le bit majoritaire de %d est res = %i \n", a, res);
		mpz_clear(tmp1);
		mpz_clear(tmp2);
		mpz_clear(tmp3);
		fmpz_poly_clear( poly_c );
	}
}
void test_bit_majoritaire(){
	unsigned a, res, aux;
	int i;
	a=11;
	mpz_t c0;
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	fhe_keygen(pk, sk);
	fmpz_poly_t poly_c;
	for(unsigned s=12; s< 32; s++){
		a=s;
		i=0;
		mpz_init(c0);
		fmpz_poly_init( poly_c );
		printf("--------->%i\n", a % 2);
		fhe_encrypt(c0, pk, a % 2);
		fmpz_poly_set_coeff_mpz ( poly_c , i , c0 );

		aux = a;	
		aux = aux >> 1;
		do {
			printf("--------->%i\n", aux % 2);
			fhe_encrypt(c0, pk, aux % 2);
			i=i+1;
			fmpz_poly_set_coeff_mpz ( poly_c , i , c0 );
			aux = aux >> 1;

		}while(aux != 0);
		mpz_t tmp1, tmp2, tmp3;
		mpz_init(tmp1);
		mpz_init(tmp2);
		mpz_init(tmp3);
		int k;
//première itération
		fmpz_poly_get_coeff_mpz ( tmp3 , poly_c , 0 );
//
		for(int j=1; j<= i; j++){
			fmpz_poly_get_coeff_mpz ( tmp2 , poly_c , j );
			fhe_mul(tmp3, tmp3, tmp2, pk);
		}
//
		for(int j= i; j>=0; j--){
			fmpz_poly_get_coeff_mpz ( tmp1 , poly_c , j );
			not(tmp1, tmp1, pk);
			k=j-1;
			while(k >=0){
				fmpz_poly_get_coeff_mpz ( tmp2 , poly_c , k );
				fhe_mul(tmp1, tmp1, tmp2, pk);
				k=k-1;
			}
			k= j+1;
			while(k <= i){
				fmpz_poly_get_coeff_mpz ( tmp2 , poly_c , k );
				fhe_mul(tmp1, tmp1, tmp2, pk);
				k=k+1;
			}
			or(tmp3, tmp1, tmp3, pk);
		}
		res = fhe_decrypt(tmp3, sk);
		printf("le bit majoritaire de %d est res = %i \n", a, res);
		mpz_clear(tmp1);
		mpz_clear(tmp2);
		mpz_clear(tmp3);
		fmpz_poly_clear( poly_c );
	}
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);
}
/*void test_sum_bits(){

	printf("FULLADD\n");
	mpz_t c0, c1, ci;
	mpz_t sum, carry;
	
	mpz_init(c0);
	mpz_init(c1);
	mpz_init(sum);
	mpz_init(carry);
	mpz_init(ci);	

	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);

	fhe_keygen(pk, sk);
	int m=15;
	fhe_encrypt(ci, pk, 0);
	fhe_encrypt(c0, pk, m % 2);
	printf("m % 2 = %u ",  m % 2);
	int aux = m / 2;

	do {
		fhe_encrypt(c1, pk, aux % 2);
		printf("aux % 2 = %u ",  aux % 2);
		aux = aux / 2;
		//mpz_set(ci,co) ;

		fhe_fulladd(sum, carry, c0, c1, ci, pk);   	
		printf("sum : %i \n",fhe_decrypt(sum, sk) );
	printf("carry : %i \n ", fhe_decrypt(carry, sk));
	mpz_set(c0, sum);
	mpz_set(ci, carry);

	}while(aux != 0);

*/


/*
	int m, m1, aux;
	mpz_t c0, c1, co, ci;
	mpz_init(c0);
	mpz_init(c1);
	mpz_init(co);
	mpz_init(ci);
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	fhe_keygen(pk, sk);
	m=15;
	aux = m;
	fhe_encrypt(c0, pk, m % 2);
	printf("**** %i\n",m % 2);
	fhe_encrypt(ci, pk, 0);

	aux =  aux >> 1;
	fhe_encrypt(c1, pk, aux % 2);
	printf("**** %i\n",aux % 2);
	//mpz_set(ci,co) ;
	fhe_fulladd(c0, co , c0, c0, ci, pk);


	do {
		printf("--------->%i\n", aux % 2);
		fhe_encrypt(c1, pk, aux % 2);
		aux = aux >> 1;
		mpz_set(ci,co) ;
		fhe_fulladd(c0, co , c0, c1, ci, pk);
	}while(aux != 0);
	
	m1 = fhe_decrypt(c0, sk);
	printf("///////////////----majorité-------!!!!!!   %i \n", m1);
	printf("decrypt co %i\n", fhe_decrypt(co, sk));


}
*/

void
test_encryt_decrypt()
{
	printf("ENCRYPT/DECRYPT\n");
	int m0, m1;
	mpz_t c0, c1;
	
	mpz_init(c0);
	mpz_init(c1);
	
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	
	for (int i = 0; i < KEYRUNS; i++) {
		fhe_keygen(pk, sk);
		
		for (int j = 0; j < RUNS; j++) {
			fhe_encrypt(c0, pk, 0);
			m0 = fhe_decrypt(c0, sk);
			fhe_encrypt(c1, pk, 1);
			m1 = fhe_decrypt(c1, sk);
			
			assert(m0 == 0);
			assert(m1 == 1);
			printf(".");
			fflush(stdout);
		}
		printf("\n");
	}
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);
	mpz_clear(c0);
	mpz_clear(c1);
	printf("PASSED.\n");
}


void
test_halfadd()
{
	printf("HALFADD\n");
	mpz_t c0, c1;
	mpz_t sum, carry;
	
	mpz_init(c0);
	mpz_init(c1);
	mpz_init(sum);
	mpz_init(carry);
	
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	
	for (int i = 0; i < KEYRUNS; i++) {
		fhe_keygen(pk, sk);
		
		fhe_encrypt(c0, pk, 0);
		fhe_encrypt(c1, pk, 1);
		
		ASSERT_HALFADD(c0,c0,0,0);
		ASSERT_HALFADD(c1,c0,1,0);
		ASSERT_HALFADD(c0,c1,1,0);
		ASSERT_HALFADD(c1,c1,0,1);
		printf(".");
		fflush(stdout);
	}
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);
	mpz_clear(sum);
	mpz_clear(carry);
	mpz_clear(c0);
	mpz_clear(c1);
	printf(" PASSED.\n");
}


void
test_fulladd()
{
	printf("FULLADD\n");
	mpz_t c0, c1;
	mpz_t sum, carry;
	
	mpz_init(c0);
	mpz_init(c1);
	mpz_init(sum);
	mpz_init(carry);
	
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	
	for (int i = 0; i < KEYRUNS; i++) {
		fhe_keygen(pk, sk);
		
		fhe_encrypt(c0, pk, 0);
		fhe_encrypt(c1, pk, 1);
		
		ASSERT_FULLADD(c0,c0,c0,0,0);
		ASSERT_FULLADD(c1,c0,c0,1,0);
		ASSERT_FULLADD(c0,c1,c0,1,0);
		ASSERT_FULLADD(c1,c1,c0,0,1);
		ASSERT_FULLADD(c0,c0,c1,1,0);
		ASSERT_FULLADD(c1,c0,c1,0,1);
		ASSERT_FULLADD(c0,c1,c1,0,1);
		ASSERT_FULLADD(c1,c1,c1,1,1);
		printf(".");
		fflush(stdout);
	}
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);
	mpz_clear(sum);
	mpz_clear(carry);
	mpz_clear(c0);
	mpz_clear(c1);
	printf(" PASSED.\n");
}


void
test_recrypt()
{
	printf("RECRYPT\n");
	
	mpz_t c0, c1;
	
	mpz_init(c0);
	mpz_init(c1);
	
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	
	for (int i = 0; i < KEYRUNS; i++) {
		fhe_keygen(pk, sk);
		
		for (int j = 0; j < RUNS; j++) {
			fhe_encrypt(c0, pk, 0);
			fhe_encrypt(c1, pk, 1);
			
			fhe_recrypt(c0, pk);
			assert(fhe_decrypt(c0, sk) == 0);
			
			fhe_recrypt(c1, pk);
			assert(fhe_decrypt(c1, sk) == 1);
			
			printf(".");
			fflush(stdout);
		}
		printf("\n");
	}
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);
	mpz_clear(c0);
	mpz_clear(c1);
	printf("PASSED.\n");
}


void
test_homomorphic()
{
	printf("HOMOMORPHIC (w/o recrypt)\n");
	
	int m;
	mpz_t c0, c1, temp;
	
	mpz_init(c0);
	mpz_init(c1);
	mpz_init(temp);
	
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	
	for (int i = 0; i < KEYRUNS; i++) {
		mpz_t c0, c1;
		
		mpz_init(c0);
		mpz_init(c1);
		
		fhe_pk_t pk;
		fhe_sk_t sk;
		fhe_pk_init(pk);
		fhe_sk_init(sk);
		
		fhe_keygen(pk, sk);
		fhe_encrypt(c0, pk, 0);
		printf("\nadd-chain: ");
		for (int j = 0; j < RUNS*RUNS; j++) {
			fhe_add(c0, c0, c0, pk);
			m = fhe_decrypt(c0, sk);
			printf("%i", m);
			fflush(stdout);
		}
		fhe_encrypt(c1, pk, 1);
		printf("\nmul-chain: ");
		for (int j = 0; j < RUNS*RUNS; j++) {
			fhe_mul(c1, c1, c1, pk);
			m = fhe_decrypt(c1, sk);
			printf("%i", m);
			fflush(stdout);
		}
		printf("\n");
	}
	
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);
	mpz_clear(c0);
	mpz_clear(c1);
	mpz_clear(temp);
	
	printf("PASSED.\n");
}

void
test_fully_homomorphic()
{
	printf("FULLY HOMOMORPHIC (with recrypt)\n");
	
	int m;
	mpz_t c0, c1, temp;
	
	mpz_init(c0);
	mpz_init(c1);
	mpz_init(temp);
	
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	
	for (int i = 0; i < KEYRUNS; i++) {
		mpz_t c0, c1;
		
		mpz_init(c0);
		mpz_init(c1);
		
		fhe_pk_t pk;
		fhe_sk_t sk;
		fhe_pk_init(pk);
		fhe_sk_init(sk);
		
		fhe_keygen(pk, sk);
		fhe_encrypt(c0, pk, 1);
		fhe_encrypt(c1, pk, 1);
		printf("\nadd-chain: ");
		for (int j = 0; j < RUNS*RUNS; j++) {
			fhe_add(c0, c0, c1, pk);
			fhe_recrypt(c0, pk);
			m = fhe_decrypt(c0, sk);
			printf("%i", m);
			fflush(stdout);
		}
		fhe_encrypt(c1, pk, 1);
		printf("\nmul-chain: ");
		for (int j = 0; j < RUNS*RUNS; j++) {
			fhe_mul(c1, c1, c1, pk);
			fhe_recrypt(c1, pk);
			m = fhe_decrypt(c1, sk);
			printf("%i", m);
			fflush(stdout);
		}
		printf("\n");
	}
	
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);
	mpz_clear(c0);
	mpz_clear(c1);
	mpz_clear(temp);
	
	printf("PASSED.\n");
}
