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


void test_keygen();

void test_suite()
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
	test_min_max();
	//test_insertion_sort();
	//test_oddeven_merger_sort();
	//test_bitonic_sort();
	//test_majority_bit();

	//test_keygen();
}

void test_keygen(){
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	fhe_keygen(pk, sk);
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);
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

void not(mpz_t res, mpz_t a, fhe_pk_t pk){
	mpz_t aux, c1;
	mpz_init(aux);
	mpz_init(c1);
	fhe_encrypt(c1, pk, 1);
	fhe_add(res, a, c1, pk);
	mpz_clear(aux);
	mpz_clear(c1);
}

void test_aIsGreater(mpz_t res, fmpz_poly_t polya, fmpz_poly_t polyb, fhe_pk_t pk, int nbits){

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


void min_max(mpz_t *min, mpz_t *max, fmpz_poly_t poly_c1, fmpz_poly_t poly_c2, fhe_pk_t pk, int nbits){
	
	mpz_t a_k;
	mpz_t b_k;
	mpz_t tmp;

	mpz_init(a_k);
	mpz_init(b_k);
	mpz_init(tmp);
	
	mpz_t aIsGreater;
	mpz_init(aIsGreater);
		
	int k;

	test_aIsGreater(aIsGreater,poly_c1,poly_c2,pk, nbits);

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
	a=1450; b=1030;  
	printf("a = %d et b = %d\n", a, b);
	aux1 = a ; aux2=b;
	int i = 0;

	int nbits;   // Number of bits in the binary representation of the integers

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

	}while(aux1 != 0 || aux2 !=0);


	/////////// Evaluation ////////////////////
	nbits= i +1;
	fmpz_poly_t max;
	fmpz_poly_t min;
	//mpz_t * max;
	//mpz_t * min;
	fmpz_poly_init(max);
	fmpz_poly_init(min);
	//max = malloc(sizeof(mpz_t) * nbits);
	//min = malloc(sizeof(mpz_t) * nbits);
	//for(i=0;i<nbits;i++){
	//	mpz_init(max[i]);
	//	mpz_init(min[i]);
	//}
	
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

		fmpz_poly_set_coeff_mpz(max , k , tmp) ;
		//mpz_set(max[k],tmp);

		fmpz_poly_get_coeff_mpz(a_k, poly_c1,k);	
		fmpz_poly_get_coeff_mpz(b_k, poly_c2,k);
			
		fhe_mul(b_k, b_k, aIsGreater,pk);
		not(tmp, aIsGreater,pk);
		fhe_mul(a_k, a_k, tmp,pk);
		or(tmp, a_k,b_k, pk);

		fmpz_poly_set_coeff_mpz(min , k , tmp) ;
		//mpz_set(min[k],tmp);		
			
	}

	
	//////////////// Decryption ///////////////////////////
	
	//for (k=0;k<nbits;k++){
	//	printf("The bit at the  %d  place is  max = %i \n", k,fhe_decrypt(max[k],sk));
	//} 
	
	//for (k=0;k<nbits;k++){
	//	printf("The bit at the  %d  place is  min = %i \n", k,fhe_decrypt(min[k],sk));
	//} 

	aux1= 0; aux2= 0;
	unsigned d;
	for(k=nbits-1; k>=0 ;k--){
		fmpz_poly_get_coeff_mpz(tmp, max ,k);
		d =  fhe_decrypt(tmp,sk);
		aux1= (aux1 * 2) + d;
		
	}
	printf("le max est: %d \n", aux1);
	for(k=nbits-1; k>=0 ;k--){
		fmpz_poly_get_coeff_mpz(tmp, min ,k);
		d= fhe_decrypt(tmp,sk);
		aux2= (aux2 * 2) +d;
	}
	printf("le min est: %d\n", aux2);

	//for(k=0;k<nbits;k++){
	//	mpz_clear(max[k]);
	//	mpz_clear(min[k]);
	//}
	
	fmpz_poly_clear( poly_c1 );
	fmpz_poly_clear( poly_c2 ); 
	
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);	
	mpz_clear(c0);
	mpz_clear(c1);
	mpz_clear(a_k);
	mpz_clear(b_k); 
	mpz_clear(tmp);
	fmpz_poly_init(max);
	fmpz_poly_init(min);
}


void OddEvenMerge(fmpz_poly_t * poly_nums, int n, int nbits, int lo, int r, fhe_sk_t sk , fhe_pk_t pk){

		
	mpz_t * max;
	mpz_t * min;
	max = malloc(sizeof(mpz_t) * nbits);
	min = malloc(sizeof(mpz_t) * nbits);
	for(int i=0;i<nbits;i++){
		mpz_init(max[i]);
		mpz_init(min[i]);
	}

	int m = 2*r;
	
	if(m<n){
		
		OddEvenMerge(poly_nums, n, nbits, lo, m, sk, pk); // even subsequence
		OddEvenMerge(poly_nums, n, nbits, lo+r, m, sk, pk); //odd subsequence	
		for(int i=lo+r; i+r <lo+n; i+=m){
			min_max(min, max, poly_nums[i], poly_nums[i+r], pk, nbits);		
			for(int k=0;k<nbits;k++){
				fmpz_poly_set_coeff_mpz(poly_nums[i], k, min[k]);
				fmpz_poly_set_coeff_mpz(poly_nums[i+r],k, max[k]);
			}

			
		}

	}	
	else {
	
		min_max(min, max, poly_nums[lo], poly_nums[lo+r], pk, nbits);		
		for(int k=0;k<nbits;k++){
			fmpz_poly_set_coeff_mpz(poly_nums[lo], k, min[k]);
			fmpz_poly_set_coeff_mpz(poly_nums[lo+r],k, max[k]);
		}


		
	}



}



void OddEvenMergeSort(fmpz_poly_t * poly_nums, int n, int nbits, int lo, fhe_sk_t sk, fhe_pk_t pk)
{
	if(n>1){
	
		int m = n/2;
		OddEvenMergeSort(poly_nums, m, nbits, lo, sk, pk);
		OddEvenMergeSort(poly_nums, m, nbits, lo+m, sk, pk);
		OddEvenMerge(poly_nums, n, nbits, lo, 1, sk, pk);
	
	}
}


void test_oddeven_merger_sort(){
	int n = 8; // Nombre d'entiers à trier
	int nbits = 6; // Size in number of bits
	int list [] = {19,8,32,13,5,6,10,24};  // List d'entiers à trier

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
	
	printf("\n");
	

	////////////////// Odd Even Sorting //////////////

	OddEvenMergeSort(poly_nums, n, nbits, 0, sk, pk);

	////////////// Decryption /////////////////////


		for(i=0;i<n;i++){
		printf(" \n The %d'th number is \n", i);
		for(k=0;k<nbits;k++){
			fmpz_poly_get_coeff_mpz(tmp, poly_nums[i],k);
			printf(" %i ", fhe_decrypt(tmp, sk));
		}
	
	}
	
	printf("\n");


}




void bitonicMergeUp(fmpz_poly_t *poly_nums, int  nbits, int n,  int lo , int  high, fhe_sk_t sk, fhe_pk_t pk){
	if(high==0) return;
	
	//////////// SWAP //////////

	mpz_t * max;
	mpz_t * min;
	max = malloc(sizeof(mpz_t) * nbits);
	min = malloc(sizeof(mpz_t) * nbits);
	for(int i=0;i<nbits;i++){
		mpz_init(max[i]);
		mpz_init(min[i]);
	}


	for(int i = 0; i<high;i++){
		min_max(min, max, poly_nums[i+lo], poly_nums[i+lo+high], pk, nbits);		
		for(int k=0;k<nbits;k++){
			fmpz_poly_set_coeff_mpz(poly_nums[i+lo], k, min[k]);
			fmpz_poly_set_coeff_mpz(poly_nums[i+lo+high],k, max[k]);
		}
		
		
	}

	bitonicMergeUp(poly_nums, nbits, n, lo, high/2,  sk, pk);
	bitonicMergeUp(poly_nums, nbits, n, lo+high, high/2,  sk, pk);
	for(int k=0;k<nbits;k++){
		mpz_clear(max[k]);
		mpz_clear(min[k]);
	}
	
				
}

void bitonicMergeDown(fmpz_poly_t *poly_nums, int  nbits, int n,  int lo , int  high, fhe_sk_t sk, fhe_pk_t pk){
	if(high==0) return;

	//////////// SWAP //////////////////

	mpz_t * max;
	mpz_t * min;
	max = malloc(sizeof(mpz_t) * nbits);
	min = malloc(sizeof(mpz_t) * nbits);
	for(int i=0;i<nbits;i++){
		mpz_init(max[i]);
		mpz_init(min[i]);
	}


	for(int i = 0; i<high;i++){
		min_max(min, max, poly_nums[i+lo], poly_nums[i+lo+high], pk, nbits);		
		for(int k=0;k<nbits;k++){
			fmpz_poly_set_coeff_mpz(poly_nums[i+lo+high], k, min[k]);
			fmpz_poly_set_coeff_mpz(poly_nums[i+lo],k, max[k]);
		}
			
	}

	bitonicMergeDown(poly_nums, nbits, n, lo, high/2,  sk, pk);
	bitonicMergeDown(poly_nums, nbits, n, lo+high, high/2,  sk, pk);
	for(int k=0;k<nbits;k++){
		mpz_clear(max[k]);
		mpz_clear(min[k]);
	}
	
}


void bitonicSortDown(fmpz_poly_t *poly_nums, int  nbits, int n,  int lo , int  high, fhe_sk_t sk, fhe_pk_t pk){
	if(high==1) return;

	bitonicSortUp(poly_nums, nbits, n, lo, high/2 , sk, pk);
	bitonicSortDown(poly_nums, nbits, n, lo+high/2, high/2, sk, pk );
	bitonicMergeDown(poly_nums, nbits, n, lo, high/2,  sk, pk);

}

void bitonicSortUp(fmpz_poly_t *poly_nums, int  nbits, int n,  int lo , int  high, fhe_sk_t sk, fhe_pk_t pk){
	if(high==1) return;

	bitonicSortUp(poly_nums, nbits, n, lo, high/2 , sk, pk);
	bitonicSortDown(poly_nums, nbits, n, lo+high/2, high/2,  sk, pk );
	bitonicMergeUp(poly_nums, nbits, n, lo, high/2,  sk, pk);

}




void test_bitonic_sort(){ // Sorts a bitonic array of 2^n elements
	int n = 8; // Nombre d'entiers à trier
	int nbits = 6; // Size in number of bits
	int list [] ={5, 2, 11, 7, 18, 4, 5, 9}; // {8,13,19,32,24,10,6,5};  // List d'entiers à trier

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
	
	printf("\n");
	

	////////////////// Bitonic Sorting : sort from the index 0 till n //////////////

	bitonicSortUp(poly_nums, nbits, n, 0, n, sk, pk);

	////////////// Decryption /////////////////////


		for(i=0;i<n;i++){
		printf(" \n The %d'th number is \n", i);
		for(k=0;k<nbits;k++){
			fmpz_poly_get_coeff_mpz(tmp, poly_nums[i],k);
			printf(" %i ", fhe_decrypt(tmp, sk));
		}
	
	}
	
	printf("\n");	
		
	for(k=0;k<n;k++)
		fmpz_poly_clear( poly_nums[k]);
	
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);	
	mpz_clear(c0);
	mpz_clear(tmp);
				


}

void test_insertion_sort(){   // Generate n random numbers of nbits each and sort them
	int n = 5; // Nombre d'entiers à trier
	int nbits = 6; // Size in number of bits
	int list [] = {19,8,32,13,5};  // List d'entiers à trier

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
			min_max(min, max, poly_nums[j], poly_nums[j+1],pk, nbits);
				
	
			for(k=0;k<nbits;k++){
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

void test_sum_mpz_t(fmpz_poly_t sum_result, fmpz_poly_t poly_c1, fmpz_poly_t poly_c2, fhe_pk_t pk , fhe_sk_t sk ){	
/*
	unsigned a ,b, aux1, aux2;
	a=29; b=96;   // Integers to be added
	aux1 = a ; aux2=b;
	int i = 0;
	int nbits = 7;   // Number of bits in the binary representation of the largest integer
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

	}while(aux1 != 0 || aux2 !=0);

*/	/////////// Evaluation ////////////////////
	int nbits;
	//printf("degree of the first polynomials is %d \n ",fmpz_poly_degree(poly_c1));
	//printf("degree of the second polynomial is %d \n ",fmpz_poly_degree(poly_c2));
	if(fmpz_poly_degree(poly_c1)>fmpz_poly_degree(poly_c2))	
		nbits = fmpz_poly_degree(poly_c1)+1;
	else nbits = fmpz_poly_degree(poly_c2)+1;
	
	//printf("nbits of the sum would be %d", nbits);
		
	/*mpz_t * sum_result;
	sum_result = malloc(sizeof(mpz_t) * (nbits+1));
	for(i=0;i<nbits+1;i++)
	mpz_init(sum_result[i]);*/

	mpz_t tmp;
	mpz_init(tmp);

	//printf("\n In the sum function : The input polynomials are : \n");
/*
	for (int j=fmpz_poly_degree(poly_c1);j>=0;j--){
		fmpz_poly_get_coeff_mpz(tmp,poly_c1,j);
		printf(" %i ", fhe_decrypt(tmp,sk));
	}
	printf(" AND \n");

	for (int j=fmpz_poly_degree(poly_c2);j>=0;j--){
		fmpz_poly_get_coeff_mpz(tmp,poly_c2,j);
		printf(" %i ", fhe_decrypt(tmp,sk));
	}
*/	

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
//		printf("for k = %d,  a_k is %i ", k, fhe_decrypt(a_k,sk));
//		printf(" b_k is %i  \n ", fhe_decrypt(b_k, sk));
		//printf("b's bit is %i \n", fhe_decrypt(b_k, sk));
		fhe_fulladd(tmp, c_out, a_k, b_k, c_in, pk);
	//	fmpz_poly_get_coeff_mpz(a_k, sum_result, k);
//		printf("sum at %d th loop is %i \n", k, fhe_decrypt(tmp, sk));
//		printf("Carry at %d th loop is %i \n", k, fhe_decrypt(c_out, sk));
		fmpz_poly_set_coeff_mpz(sum_result,k,tmp);
		mpz_set(c_in, c_out);
					
	}
	
	fmpz_poly_set_coeff_mpz(sum_result, nbits, c_out);
	
	//////////////// Decryption ///////////////////////////
	//printf("\n The sum is \n");

	//for (k=fmpz_poly_degree(sum_result);k>=0;k--){
	//	fmpz_poly_get_coeff_mpz(tmp, sum_result, k);
		//printf(" %i ",fhe_decrypt(tmp,sk));
	//} 
	//printf("\n Exiting Sum function \n");
	//for(k=0;k<nbits+1;k++)
	//	mpz_clear(sum_result[k]);

	//free(sum_result);
	//fmpz_poly_clear( poly_c1 );
	//fmpz_poly_clear( poly_c2 ); 
	
	
	mpz_clear(c_in);
	mpz_clear(c_out);	
	mpz_clear(tmp);
	mpz_clear(a_k);
	mpz_clear(b_k); 
}

void test_majority_bit(){
	
	unsigned a ,aux1, aux2;
	a=35;  // 
	aux1= a; 
	int i = 0;
	int nbits=6;
	aux2=nbits;
	mpz_t c0;
	mpz_t c1;
	
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
	fhe_encrypt(c1, pk, nbits % 2);
	
	fmpz_poly_set_coeff_mpz( poly_c1 , i , c0 );
	fmpz_poly_set_coeff_mpz(poly_c2, i, c1);		
	aux1 = aux1 >> 1;
	aux2 = aux2 >> 1;
	do {
		//printf("--------->%i\n", aux % 2);
		fhe_encrypt(c0, pk, aux1 % 2);
		fhe_encrypt(c1, pk, aux2 %2);
		i++;
		fmpz_poly_set_coeff_mpz (poly_c1 , i , c0 );
		fmpz_poly_set_coeff_mpz (poly_c2, i, c1);
		aux1 = aux1 >> 1;
		aux2 = aux2 >> 1;
	
	}while(aux1 != 0 || aux2!=0 );

	mpz_t tmp;
	mpz_init(tmp);	

	////////////////// Evaluation ///////////////////////

	printf("Check After encryption \n");
	for (int k=nbits-1;k>=0;k--){
		fmpz_poly_get_coeff_mpz(tmp, poly_c1,k);
		printf(" %i ", fhe_decrypt(tmp,sk));
	}

	printf("\n");		
	
	fmpz_poly_t tmp1_poly;
	fmpz_poly_init(tmp1_poly);
	fmpz_poly_t tmp2_poly;
	fmpz_poly_init(tmp2_poly);


	mpz_t tmp1 ;
	mpz_init(tmp1);
	mpz_t tmp2 ;
	mpz_init(tmp2);
/*	

	fhe_encrypt(tmp1, pk, 0);
	fmpz_poly_set_coeff_mpz(tmp1_poly,0,tmp1);
	//fmpz_poly_set_coeff_mpz(tmp2_poly,0,tmp1);
	//fhe_encrypt(tmp1,pk,1);
	//fmpz_poly_set_coeff_mpz(tmp2_poly,1,tmp1);
	//test_sum_mpz_t(tmp1_poly, tmp2_poly, tmp1_poly,pk,sk);
	
	
	for(long k = 0 ; k<nbits; k++){		
		fmpz_poly_get_coeff_mpz(tmp2, poly_c1, k);
		fmpz_poly_set_coeff_mpz(tmp2_poly,0,tmp2);
		
		printf("\n the sum of %dth bit \n", k);
		for (int j=fmpz_poly_degree(tmp2_poly);j>=0;j--){
			fmpz_poly_get_coeff_mpz(tmp1,tmp2_poly,j);
			printf(" %i ", fhe_decrypt(tmp1,sk));
		}

		printf("\n AND the temporary sum \n");

		for (int j=fmpz_poly_degree(tmp1_poly);j>=0;j--){
			fmpz_poly_get_coeff_mpz(tmp1,tmp1_poly,j);
			printf(" %i ", fhe_decrypt(tmp1,sk));
		}		

		printf("\n IS \n ");
	
		test_sum_mpz_t(tmp1_poly, tmp1_poly, tmp2_poly,pk, sk);

		for (int j=fmpz_poly_degree(tmp1_poly);j>=0;j--){
			fmpz_poly_get_coeff_mpz(tmp1,tmp1_poly,j);
			printf(" %i ", fhe_decrypt(tmp1,sk));
		}

				
	}
	
	printf(" \n Printing the sum of n bits \n");

	for(int k=fmpz_poly_degree(tmp1_poly);k>=0;k--){

		fmpz_poly_get_coeff_mpz(tmp1,tmp1_poly,k);
		printf(" %i ",fhe_decrypt(tmp1,sk));
	}
	

	//// Calling isGreaterfunction //////////

	test_aIsGreater(tmp1,tmp1_poly, poly_c2, pk, fmpz_poly_degree(tmp1_poly)+1);

	printf("The majority bit in %d is %i", a, fhe_decrypt(tmp1, sk) ) ;

*/	
	
	printf("\n");
	fmpz_poly_clear(poly_c1);
	fmpz_poly_clear(poly_c2);
	mpz_clear(c1);
	fmpz_poly_clear(tmp1_poly);
	fmpz_poly_clear(tmp2_poly);
	mpz_clear(tmp1);
	mpz_clear(tmp2);
	fhe_pk_clear(pk);
	mpz_clear(c0);
	fhe_sk_clear(sk);
	mpz_clear(tmp);
}








void test_sum_integers(){	

	unsigned a ,b, aux1, aux2;
	a=29; b=96;   // Integers to be added
	aux1 = a ; aux2=b;
	int i = 0;
	int nbits = 7;   // Number of bits in the binary representation of the largest integer
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

	}while(aux1 != 0 || aux2 !=0);

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
		printf("for k = %d,  a_k is %i \n", k, fhe_decrypt(a_k,sk));
		printf("for k = %d, b_k is %i  \n ", k, fhe_decrypt(b_k, sk));
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
