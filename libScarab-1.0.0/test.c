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
	//test_xor_bits();
	//test_sum_bits();
	//test_bit_majoritaire();
	//test_sum_integers();
	//test_min_max();
	test_insertion_sort();
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

void test_aIsGreater(mpz_t res, fmpz_poly_t polya, fmpz_poly_t polyb, fhe_pk_t pk){

	long nbits = fmpz_poly_degree(polya)+1;
	mpz_t a_k;
	mpz_t b_k;
	mpz_t tmp1;
	mpz_t anot;
	mpz_t bnot;
	mpz_t tmp2;
	mpz_t count;

	mpz_init(a_k);
	mpz_init(b_k);
	mpz_init(tmp1);
	mpz_init(tmp2);
	mpz_init(anot);
	mpz_init(bnot);
	mpz_init(count);	

	mpz_t aIsGreater;
	mpz_init(aIsGreater);

	fhe_encrypt(aIsGreater, pk, 0);	
	fhe_encrypt(count, pk, 1);

	for(long k=nbits-1;k>=0;k--){

		fmpz_poly_get_coeff_mpz(a_k, polya,k);	
		fmpz_poly_get_coeff_mpz(b_k, polyb,k);
		
		not(anot, a_k, pk);
		not(bnot, b_k, pk);	

		fhe_mul(tmp1, a_k, b_k, pk);
		fhe_mul(tmp2, anot, bnot,pk); 
		or(tmp1, tmp2, tmp1, pk); //Ek

		fhe_mul(tmp2, a_k, bnot, pk);
		fhe_mul(tmp2, tmp2, count, pk);		
		
		fhe_mul(count, count, tmp1,pk);		
		or(aIsGreater,aIsGreater, tmp2,pk);		
	}
	mpz_set(res,aIsGreater);
	mpz_clear(anot);
	mpz_clear(bnot);
	mpz_clear(count);
	mpz_clear(a_k);
	mpz_clear(b_k);
	mpz_clear(tmp1);
	mpz_clear(tmp2);
	mpz_clear(aIsGreater);
}


void min_max(mpz_t *min, mpz_t *max, fmpz_poly_t poly_c1, fmpz_poly_t poly_c2, fhe_pk_t pk){
	

	int nbits = fmpz_poly_degree(poly_c1)+1; 
	
	mpz_t a_k;
	mpz_t b_k;
	mpz_t tmp;

	mpz_init(a_k);
	mpz_init(b_k);
	mpz_init(tmp);
	
	mpz_t aIsGreater;
	mpz_init(aIsGreater);
		
	int k;

	test_aIsGreater(aIsGreater,poly_c1,poly_c2,pk);

	//printf("Is a greater than b ? Ans : %i  \n", fhe_decrypt(aIsGreater,sk));
	
	for(k=0;k<nbits;k++){
		
		fmpz_poly_get_coeff_mpz(a_k, poly_c1,k);	
		fmpz_poly_get_coeff_mpz(b_k, poly_c2,k);
		
		fhe_mul(a_k, a_k, aIsGreater,pk);
		not(tmp, aIsGreater,pk);
		fhe_mul(b_k, b_k, tmp,pk);
		or(tmp, a_k,b_k, pk);	

		mpz_set(max[k],tmp);

		fmpz_poly_get_coeff_mpz(a_k, poly_c1,k);	
		fmpz_poly_get_coeff_mpz(b_k, poly_c2,k);
			

		fhe_mul(b_k, b_k, aIsGreater,pk);
		not(tmp, aIsGreater,pk);
		fhe_mul(a_k, a_k, tmp,pk);
		or(tmp, a_k,b_k, pk);

		mpz_set(min[k],tmp);		
			
	}

	
	//////////////// Decryption ///////////////////////////

	printf("Inside min-max method \n");
	
	mpz_clear(a_k);
	mpz_clear(b_k); 
	mpz_clear(tmp);
	mpz_clear(aIsGreater);

}




void test_min_max(){
	
	unsigned a ,b, aux1, aux2;
	a=80; b=120;   // Integers to be compared suppposed to be of the same size
	aux1 = a ; aux2=b;
	int i = 0;
	int nbits = 7;   // Number of bits in the binary representation of the integers
	mpz_t c0, c1;
	
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	fhe_keygen(pk, sk);
	
	fmpz_poly_t poly_c1;
	fmpz_poly_t poly_c2;	
	
	mpz_init(c0);
	mpz_init(c1);

	fmpz_poly_init(poly_c1);
	fmpz_poly_init(poly_c2);

	////// Encryption of the bit sequences ////////////////
	fhe_encrypt(c0, pk, a % 2);
	fhe_encrypt(c1, pk, b % 2);

	fmpz_poly_set_coeff_mpz( poly_c1 , i , c0 );
	fmpz_poly_set_coeff_mpz( poly_c2, i, c1 );

	aux1 = aux1 >> 1;
	
	aux2 = aux2 >> 1;
	
	do {
		//printf("--------->%i\n", aux % 2);
		fhe_encrypt(c0, pk, aux1 % 2);
		fhe_encrypt(c1, pk, aux2 %2);
		i++;
		fmpz_poly_set_coeff_mpz ( poly_c1 , i , c0 );
		fmpz_poly_set_coeff_mpz (poly_c2, i, c1);
		aux1 = aux1 >> 1;
		aux2 = aux2 >> 1;

	}while(aux1 != 0 && aux2 !=0);


	/////////// Evaluation ////////////////////
	
	mpz_t * max;
	mpz_t * min;
	max = malloc(sizeof(mpz_t) * nbits);
	min = malloc(sizeof(mpz_t) * nbits);
	for(i=0;i<nbits;i++){
		mpz_init(max[i]);
		mpz_init(min[i]);
	}
	
	mpz_t a_k;
	mpz_t b_k;
	mpz_t tmp;

	mpz_init(a_k);
	mpz_init(b_k);
	mpz_init(tmp);
	
	mpz_t aIsGreater;
	mpz_init(aIsGreater);

	fhe_encrypt(aIsGreater, pk, 0);
		
	int k;

	for(k=0;k<nbits;k++){

		fmpz_poly_get_coeff_mpz(a_k, poly_c1,k);	
		fmpz_poly_get_coeff_mpz(b_k, poly_c2,k);
		
		not(b_k, b_k, pk);	
		
		fhe_mul(tmp, a_k, b_k, pk);
		or(aIsGreater,aIsGreater,tmp , pk);			
	}

	printf("Is a greater than b ? Ans : %i  \n", fhe_decrypt(aIsGreater,sk));
	
	for(k=0;k<nbits;k++){
		
		fmpz_poly_get_coeff_mpz(a_k, poly_c1,k);	
		fmpz_poly_get_coeff_mpz(b_k, poly_c2,k);
		
		fhe_mul(a_k, a_k, aIsGreater,pk);
		not(tmp, aIsGreater,pk);
		fhe_mul(b_k, b_k, tmp,pk);
		or(tmp, a_k,b_k, pk);	

		mpz_set(max[k],tmp);

		fmpz_poly_get_coeff_mpz(a_k, poly_c1,k);	
		fmpz_poly_get_coeff_mpz(b_k, poly_c2,k);
			

		fhe_mul(b_k, b_k, aIsGreater,pk);
		not(tmp, aIsGreater,pk);
		fhe_mul(a_k, a_k, tmp,pk);
		or(tmp, a_k,b_k, pk);

		mpz_set(min[k],tmp);		
			
	}

	
	//////////////// Decryption ///////////////////////////
	
	for (k=0;k<nbits;k++){
		printf("The bit at the  %d  place is  max = %i \n", k,fhe_decrypt(max[k],sk));
	} 
	
	for (k=0;k<nbits;k++){
		printf("The bit at the  %d  place is  min = %i \n", k,fhe_decrypt(min[k],sk));
	} 


	for(k=0;k<nbits;k++){
		mpz_clear(max[k]);
		mpz_clear(min[k]);
	}
	
	fmpz_poly_clear( poly_c1 );
	fmpz_poly_clear( poly_c2 ); 
	
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);	
	mpz_clear(c0);
	mpz_clear(c1);
	mpz_clear(a_k);
	mpz_clear(b_k); 
	mpz_clear(tmp);

}


void test_insertion_sort(){   // Generate n random numbers of nbits each and sort them
	int n = 4; // Nombre d'entiers à trier
	int nbits = 4; // Nombre de bit de chaque entier
	int list [] = {8,13,11,9};  // List d'entiers à trier

	fmpz_poly_t * poly_nums;
	poly_nums = malloc(sizeof(fmpz_poly_t) * n);
	
	int aux ; 	
	int a ;

	int i;

	mpz_t c0;
	
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	fhe_keygen(pk, sk);
	
	mpz_init(c0);


	////////////// Encryption of the bit sequennces //////////
	for(int k = 0; k < n ; k++){
		i=0;
		fmpz_poly_init(poly_nums[k]);
		a= list[k];
		aux=a;
		fhe_encrypt(c0, pk, a % 2);
		fmpz_poly_set_coeff_mpz( poly_nums[k] , i , c0 );
		aux = aux >> 1;
		do {
		//printf("--------->%i\n", aux % 2);
			fhe_encrypt(c0, pk, aux % 2);
			i++;
			fmpz_poly_set_coeff_mpz ( poly_nums[k] , i , c0 );
			aux = aux >> 1;
		}while(aux!= 0);

	}

	mpz_t tmp;	
	mpz_init(tmp);
	
	int k;
	for(i=0;i<n;i++){
		printf(" \n The %d'th number is \n", i);
		for(k=0;k<nbits;k++){
			fmpz_poly_get_coeff_mpz(tmp, poly_nums[i],k);
			printf(" %i ", fhe_decrypt(tmp, sk));
	
	}
	
	}
	


	/////////// Evaluation ////////////////////
		
	
	mpz_t * max;
	mpz_t * min;
	max = malloc(sizeof(mpz_t) * nbits);
	min = malloc(sizeof(mpz_t) * nbits);
	for(i=0;i<nbits;i++){
		mpz_init(max[i]);
		mpz_init(min[i]);
	}
	

	i=n-2;
	
	while(i>=0){
		for(int j=0;j<=i;j++){
			min_max(min, max, poly_nums[j], poly_nums[j+1],pk);
			for(k=0;k<n;k++){
				fmpz_poly_set_coeff_mpz(poly_nums[j], k, min[k]);
				fmpz_poly_set_coeff_mpz(poly_nums[j+1],k, max[k]);
			}
		
	}
	i--;
	
	}


	
	//////////////// Decryption ///////////////////////////
	
	for (i=0;i<n;i++){
		printf(" \n the %dth integer is : \n ", i);	
		for (k=0;k<nbits;k++){
			fmpz_poly_get_coeff_mpz(tmp, poly_nums[i],k);
			printf(" %i ", fhe_decrypt(tmp,sk));
		} 
		
	}
	
	printf("\n");

	for(k=0;k<nbits;k++){
		mpz_clear(max[k]);
		mpz_clear(min[k]);
	}
	
	for(k=0;k<n;k++)
		fmpz_poly_clear( poly_nums[k]);
	
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);	
	mpz_clear(c0);
	mpz_clear(tmp);
	
}

void test_sum_integers(){	

	unsigned a ,b, aux1, aux2;
	a=1692; b=1668;   // Integers to be added
	aux1 = a ; aux2=b;
	int i = 0;
	int nbits = 11;   // Number of bits in the binary representation of the integers
	mpz_t c0, c1;
	
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	fhe_keygen(pk, sk);
	
	fmpz_poly_t poly_c1;
	fmpz_poly_t poly_c2;	
	
	mpz_init(c0);
	mpz_init(c1);

	fmpz_poly_init(poly_c1);
	fmpz_poly_init(poly_c2);

	////// Encryption of the bit sequences ////////////////
	fhe_encrypt(c0, pk, a % 2);
	fhe_encrypt(c1, pk, b % 2);

	fmpz_poly_set_coeff_mpz( poly_c1 , i , c0 );
	fmpz_poly_set_coeff_mpz( poly_c2, i, c1 );
	
	aux1 = aux1 >> 1;
	
	aux2 = aux2 >> 1;
	
	do {
		//printf("--------->%i\n", aux % 2);
		fhe_encrypt(c0, pk, aux1 % 2);
		fhe_encrypt(c1, pk, aux2 %2);
		i++;
		fmpz_poly_set_coeff_mpz ( poly_c1 , i , c0 );
		fmpz_poly_set_coeff_mpz (poly_c2, i, c1);
		aux1 = aux1 >> 1;
		aux2 = aux2 >> 1;

	}while(aux1 != 0 && aux2 !=0);

	/////////// Evaluation ////////////////////
	
	mpz_t * sum_result;
	sum_result = malloc(sizeof(mpz_t) * (nbits+1));
	for(i=0;i<nbits+1;i++)
	mpz_init(sum_result[i]);
	mpz_t a_k;
	mpz_t b_k;
	mpz_t c_in;
	mpz_t c_out;

	mpz_init(a_k);
	mpz_init(b_k);
	mpz_init(c_in);	
	mpz_init(c_out);
	fhe_encrypt(c_in, pk, 0);
	fhe_encrypt(c_out, pk, 0);

	int k;
	for(k=0;k<nbits;k++){

		fmpz_poly_get_coeff_mpz(a_k, poly_c1,k);	
		fmpz_poly_get_coeff_mpz(b_k, poly_c2,k);
		//printf("b's bit is %i \n", fhe_decrypt(b_k, sk));
		fhe_fulladd(sum_result[k], c_out, a_k, b_k, c_in, pk);
		//printf("sum at %d th loop is %i \n", k, fhe_decrypt(sum_result[k], sk));
		//printf("Carry at %d th loop is %i \n", k, fhe_decrypt(c_out, sk));
		mpz_set(c_in, c_out);
					
	}
	
	mpz_set(sum_result[nbits],c_out);
	
	
	//////////////// Decryption ///////////////////////////
	
	for (k=0;k<nbits+1;k++){
		printf("The bit at the  %d  place is  sum_res = %i \n", k,fhe_decrypt(sum_result[k],sk));
	} 

	for(k=0;k<nbits+1;k++)
		mpz_clear(sum_result[k]);

	//free(sum_result);
	fmpz_poly_clear( poly_c1 );
	fmpz_poly_clear( poly_c2 ); 
	
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);
	mpz_clear(c_in);
	mpz_clear(c_out);	
	mpz_clear(c0);
	mpz_clear(c1);
	mpz_clear(a_k);
	mpz_clear(b_k); 
}


void test_xor_bits(){
	int m, aux;
	mpz_t c0, c1;
	mpz_init(c0);
	mpz_init(c1);
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	fhe_keygen(pk, sk);
	m= 13;	
	for(int i=13; i < 27; i++){
		m=i;
		printf("--------->%i\n", m % 2);
		fhe_encrypt(c0, pk, m % 2);
		aux = m;	
		aux = aux >> 1;
		do {
			printf("--------->%i\n", aux % 2);
			fhe_encrypt(c1, pk, aux % 2);
			aux = aux >> 1;
			fhe_add(c0, c0, c1, pk);
		}while(aux != 0);
////
		aux = fhe_decrypt(c0, sk);
		printf("///////////////-----------!!!!!!  xor bits de %i est %i \n", m, aux);
	}
	mpz_clear(c0);
	mpz_clear(c1);
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


/*
void test_sum_bits(){

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
