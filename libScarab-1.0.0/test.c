/*
 *  test.c
 *  integer-fhe
 *
 *  Created by Henning Perl on 17.12.10.
 *
 */

#include "test.h"
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>

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

void mesure_sum_integers();
void mesure_encrypt_decrypt();
//void mesure_decrypt();
void measure_test_keygen(){
	for(int i=0; i< 40; i++) test_keygen();
}

void bitonicSortUp(fmpz_poly_t *poly_nums, int  nbits, int n,  int lo , int  high, fhe_sk_t sk, fhe_pk_t pk);


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
	//mesure_sum_integers();
	//test_min_max();
	//test_insertion_sort();
	//test_oddeven_merger_sort();
	//test_bitonic_sort();
	//test_majority_bit();

	//mesure_test_matrix_prod();
	//test_sum_unbounded();
	//measure_test_keygen();
	//test_keygen();
	mesure_encrypt_decrypt();
	//essai();
}

void essai(){
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	fhe_keygen(pk, sk);
	mpz_t c;
	mpz_init(c);
	mpz_set_ui(c, 1);
	int m= fhe_decrypt(c, sk);
	printf("m= %d\n",m);
	mpz_clear(c);
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);
}
void test_keygen(){
	struct timeval start, end;	
    	long mtime, seconds, useconds;    
	//clock_t  START_eval;
	//double T_Elapsed4;
	//START_eval = clock();
 	gettimeofday(&start, NULL);    
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	fhe_keygen(pk, sk);
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);
	//T_Elapsed4 = (double) (clock () - START_eval);
	//printf(" Evaluation took %f clock/sec \n ", T_Elapsed4);     	
	gettimeofday(&end, NULL);
	seconds  = end.tv_sec  - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
	printf("Elapsed time in KeyGen : %ld ms \n", mtime );
}

void mesure_encrypt_decrypt(){
	fhe_pk_t pk;
	fhe_sk_t sk;	
	mpz_t c, aux1, aux2;
	int m;
	mpz_init(c);
	mpz_init(aux1);
	mpz_init(aux2);
	mpz_set_ui(aux2, 0);
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	struct timeval start, end;	
    	long mtime, seconds, useconds;  
	gettimeofday(&start, NULL);
	fhe_keygen(pk, sk);
	gettimeofday(&end, NULL);
	seconds  = end.tv_sec  - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
	//gmp_printf("p: %Zd \n", pk->p);
	//gmp_printf("alpha: %Zd \n", pk->alpha);
	//gmp_printf("B: %Zd \n", sk->B);

	printf("N: %d , log(p): %d , Mu: %d , Nu: %d , keygen: %ld ",N,mpz_sizeinbase( pk->p, 2), MU, LOG_NU, mtime );
	
	
	gettimeofday(&start, NULL);  
	for(int i=0; i < 100; i++){
		fhe_encrypt(c, pk, 1);
	}
	gettimeofday(&end, NULL);
	seconds  = end.tv_sec  - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
	printf("encrypt: %ld ms ", mtime );
	

	gettimeofday(&start, NULL); 
	for(int i=0; i < 100; i++){   
		fhe_fulladd(c, aux1,c,c,aux2, pk);
	}
	gettimeofday(&end, NULL);
	seconds  = end.tv_sec  - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
	printf("add : %ld ms ", mtime );


	gettimeofday(&start, NULL);  
	for(int i=0; i < 100; i++){  
		fhe_mul(c,c, c, pk);
	}
	gettimeofday(&end, NULL);
	seconds  = end.tv_sec  - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
	printf("mul : %ld ms ", mtime );
	
	gettimeofday(&start, NULL);  
	for(int i=0; i < 100; i++){
		fhe_recrypt(c, pk);
	}
	gettimeofday(&end, NULL);
	seconds  = end.tv_sec  - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
	printf("recrypt : %ld ms ", mtime );
	
	gettimeofday(&start, NULL);  
	for(int i=0; i < 100; i++){  
		fhe_decrypt(c, sk);
	}
	//printf("m : %d ",m);
	gettimeofday(&end, NULL);
	seconds  = end.tv_sec  - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
	printf("decrypt : %ld ms ", mtime );
	
	
	gettimeofday(&start, NULL);
	test_sum_integers(2000 , 49, pk, sk);
	gettimeofday(&end, NULL);
	seconds  = end.tv_sec  - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
	printf("sum: %ld ms \n", mtime );
	mpz_clear(c); mpz_clear(aux1); mpz_clear(aux2);
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);
}

//void mesure_decrypt(){

//}

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

	//printf("Inside min-max method \n");
	
	mpz_clear(a_k);
	mpz_clear(b_k); 
	mpz_clear(tmp);
	mpz_clear(aIsGreater);

}


void test_min_max(){
	
	clock_t START_init = clock();

	struct timeval start, end;
	
    	long mtime, seconds, useconds;    

   	 gettimeofday(&start, NULL);
    	
   
 	
	////////////////  Initialization ////////////////

	unsigned a ,b, aux1, aux2;

	a=2; b=5;  

	printf("a = %d et b = %d\n", a, b);
	aux1 = a ; aux2=b;
	int i = 0;
//	unsigned a =30050183, b= 504195648;
//	unsigned aux1, aux2;


	int nbits;   // Number of bits in the binary representation of the integers

	mpz_t c0, c1;
	fmpz_poly_t poly_c1;
	fmpz_poly_t poly_c2;	
	
	mpz_init(c0);
	mpz_init(c1);

	fmpz_poly_init(poly_c1);
	fmpz_poly_init(poly_c2);

		
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	
	double T_Elapsed1 = (double) ( clock () - START_init ); 
	printf(" Initialization of the variables etc took %f clocks / sec \n ", T_Elapsed1);

	gettimeofday(&end, NULL);
	seconds  = end.tv_sec  - start.tv_sec;
  	useconds = end.tv_usec - start.tv_usec;
   	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;

   	 printf("Elapsed time in Init : %ld milliseconds\n", mtime);
  	
	////////////////////// Initialization Ends ////////////////////
	

	//////////////////////// Key Generation /////////
	
	clock_t  START_keygen = clock();
	gettimeofday(&start, NULL);


	fhe_keygen(pk, sk);

	double T_Elapsed2 = (double) (clock () - START_keygen);	
	printf(" KeyGen took %f clocks/sec \n", T_Elapsed2);
	
	gettimeofday(&end, NULL);
	seconds  = end.tv_sec  - start.tv_sec;
  	useconds = end.tv_usec - start.tv_usec;
   	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;

   	printf("Elapsed time in KeyGen : %ld milliseconds\n", mtime);
  		
	////////////////////// Key Generation Ends /////////////
	
	////// Encryption of the bit sequences ////////////////

	clock_t  START_enc = clock();
	gettimeofday(&start, NULL);

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

	nbits=i+1;

	printf("Maximum number of bits is %d", nbits);
	
	double T_Elapsed3 = (double) (clock () - START_enc);
	printf(" Encryption took %f clocks/sec \n ", T_Elapsed3);	
	gettimeofday(&end, NULL);
	seconds  = end.tv_sec  - start.tv_sec;
  	useconds = end.tv_usec - start.tv_usec;
   	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;

   	printf("Elapsed time in Encryption  : %ld milliseconds\n", mtime);
  	/////////////////// Encryption Ends /////////////



	/////////// Evaluation ////////////////////
	clock_t  START_eval = clock();
	gettimeofday(&start, NULL);

	mpz_t * max;
	mpz_t * min;
	max = malloc(sizeof(mpz_t) * nbits);
	min = malloc(sizeof(mpz_t) * nbits);
	for(i=0;i<nbits;i++){
		mpz_init(max[i]);
		mpz_init(min[i]);
	}


	/////////// Evaluation ////////////////////
	//nbits= i +1;
	//fmpz_poly_t max;
	//fmpz_poly_t min;
	//mpz_t * max;
	//mpz_t * min;
	//fmpz_poly_init(max);
	//fmpz_poly_init(min);
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

	
	min_max(min, max, poly_c1, poly_c2, pk, nbits);
	
	double T_Elapsed4 = (double) (clock () - START_eval);
	printf(" Evaluation took %f clock/sec \n ", T_Elapsed4);
					
	gettimeofday(&end, NULL);
	seconds  = end.tv_sec  - start.tv_sec;
  	useconds = end.tv_usec - start.tv_usec;
   	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;


   	printf("Elapsed time in Evaluation : %ld milliseconds\n", mtime);
  	


	//////////////// Evaluation Ends ////////////////

/*
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
			
	}*/



	///////////////////// Decryption /////////////////
	

	clock_t  START_dec = clock();
	gettimeofday(&start, NULL);

	aux1= 0; aux2= 0;

	unsigned d; int k;
	for(k=nbits-1; k>=0 ;k--){
		d =  fhe_decrypt(max[k],sk);
		aux1= (aux1 * 2) + d;
		
	}

	printf("le max est: %d \n", aux1);

	for(k=nbits-1;k>=0 ;k--){
		d= fhe_decrypt(min[k],sk);
		aux2= (aux2 * 2) +d;
	}
	printf("le min est: %d\n", aux2);
	
	double T_Elapsed5 = (double) (clock () - START_dec);
	printf(" Decryption took  %f clocks/sec \n ", T_Elapsed5);
	gettimeofday(&end, NULL);
	seconds  = end.tv_sec  - start.tv_sec;
  	useconds = end.tv_usec - start.tv_usec;
   	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;


   	printf("Elapsed time in Decryption : %ld milliseconds\n", mtime);
  	
	//////////////////////// Decryption Ends /////////////
	
	
	for(k=0;k<nbits;k++){
		mpz_clear(max[k]);
		mpz_clear(min[k]);
	}
	

	free(max);
	free(min);
	fmpz_poly_clear( poly_c1 );
	fmpz_poly_clear( poly_c2 ); 
	
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);	
	mpz_clear(c0);
	mpz_clear(c1);
	mpz_clear(a_k);
	mpz_clear(b_k); 
	mpz_clear(tmp);

	mpz_clear(aIsGreater);


}


void OddEvenMerge(fmpz_poly_t * poly_nums, int n, int nbits, int lo, int high, int r, fhe_pk_t pk){

	//printf("Inside OddEven Merge");		
	mpz_t * max;
	mpz_t * min;
	max = malloc(sizeof(mpz_t) * nbits);
	min = malloc(sizeof(mpz_t) * nbits);
	for(int i=0;i<nbits;i++){
		mpz_init(max[i]);
		mpz_init(min[i]);
	}

	
	int m = 2*r;
	
	if(m<high){
		OddEvenMerge(poly_nums, n, nbits, lo, high, m, pk); // even subsequence
		OddEvenMerge(poly_nums, n, nbits, lo+r, high, m, pk); //odd subsequence	
	for(int i=lo+r; i+r<lo+high; i+=m){
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

	for(int k=0;k<nbits;k++){
		mpz_clear(max[k]);
		mpz_clear(min[k]);
	}
	free(max);
	free(min);
	
	
}



void OddEvenMergeSort(fmpz_poly_t * poly_nums, int n, int nbits, int lo, int high , fhe_pk_t pk, fhe_sk_t sk)
{
	if(high>1){
		int m = high/2;	
		OddEvenMergeSort(poly_nums, n, nbits, lo, m, pk, sk);
		OddEvenMergeSort(poly_nums, n, nbits, lo+m, m, pk, sk);		      	      OddEvenMerge(poly_nums, n, nbits, lo, high , 1,  pk);
		
	}
	
	
}


void test_oddeven_merger_sort(){
	int n = 32;// Nombre d'entiers à trier
	int nbits = 6;// Size in number of bits
	int *list = malloc(sizeof(int)*n); // List d'entiers à trier
	//int list[]= {8,7,13,2};
	fmpz_poly_t * poly_nums;
	poly_nums = malloc(sizeof(fmpz_poly_t) * n);
	
	int aux ; 	
	int a ;

	int i;

	mpz_t c0;
	mpz_t tmp;	
	mpz_init(tmp);

	int d ;
	
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	fhe_keygen(pk, sk);
	
	mpz_init(c0);
	struct timeval start, end;
	long mtime, seconds, useconds;    
	
	for(int obs=0;obs<5;obs++){
	////////////// Encryption of the bit sequennces //////////
		srand(time(NULL));
		int mod = (int)pow(2, nbits);
		for(int t=0;t<n;t++){
			list[t]= rand()%mod;
			printf("%d ",list[t]);	
		}
	printf("\n");
	nbits=0;
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
		if(i+1>nbits) nbits=i+1;
	}

	
	////////////////// Odd Even Sorting //////////////
	clock_t  START_eval = clock();
	gettimeofday(&start, NULL);

	OddEvenMergeSort(poly_nums, n, nbits, 0, n, pk,sk);
	double T_Elapsed4 = (double) (clock () - START_eval);
		//T_Elapsed4/=CLOCKS_PER_SEC;
		//printf(" Evaluation took %f clock/sec \n ", T_Elapsed4);
					
	gettimeofday(&end, NULL);
	seconds  = end.tv_sec  - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
   	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
		
		
  	printf("%d \t %d \t %f \t %ld \n",n, nbits, T_Elapsed4, mtime);	


	
	////////////// Decryption /////////////////////
	printf("After OddEven sort \n");
	for (i=0;i<n;i++){
		aux = 0;
		d=0;
		for (int k=nbits-1;k>=0;k--){
			fmpz_poly_get_coeff_mpz(tmp, poly_nums[i],k);
			d = fhe_decrypt(tmp, sk);
			aux = (aux*2)+d;
		} 
		printf("%d ", aux);
		
	}
	printf("\n");


}

	for(int k=0;k<n;k++)
		fmpz_poly_clear( poly_nums[k]);
	free(poly_nums);
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);	
	mpz_clear(c0);
	mpz_clear(tmp);
	free(list);			



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
	free(max);
	free(min);
				
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
	free(max);
	free(min);	
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
	int n = 32;// Nombre d'entiers à trier
	int nbits = 8 ;// Size in number of bits
	int *list = malloc(sizeof(int)*n); // List d'entiers à trier
	
	fmpz_poly_t * poly_nums;
	poly_nums = malloc(sizeof(fmpz_poly_t) * n);
	
	int aux ; 	
	int a ;

	int i;

	mpz_t c0;
	mpz_t tmp;	
	mpz_init(tmp);

	int d ;
		
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	fhe_keygen(pk, sk);
	
	mpz_init(c0);
	struct timeval start, end;
	long mtime, seconds, useconds;    
	for(int obs=0;obs<3;obs++){
	////////////// Encryption of the bit sequennces //////////
		srand(time(NULL));
		int mod = (int)pow(2, nbits);
		printf("The array is \n");

		for(int t=0;t<n;t++){
			list[t]= rand()%mod;
			printf("%d ", list[t]);
			
		}
		printf("\n");
		nbits = 0;
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
			if(i+1>nbits) nbits=i+1;
		}


	////////////////// Bitonic Sorting : sort from the index 0 till n //////////////
		clock_t  START_eval = clock();
		gettimeofday(&start, NULL);

		bitonicSortUp(poly_nums, nbits, n, 0, n, sk, pk);
		
		double T_Elapsed4 = (double) (clock () - START_eval);
		//T_Elapsed4/=CLOCKS_PER_SEC;
		//printf(" Evaluation took %f clock/sec \n ", T_Elapsed4);
					
		gettimeofday(&end, NULL);
		seconds  = end.tv_sec  - start.tv_sec;
		useconds = end.tv_usec - start.tv_usec;
   		mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
		
		
  	 	printf("%d \t %d \t %f \t %ld \n",n, nbits, T_Elapsed4, mtime);	
		printf("After Bitonic sort \n");
		for (i=0;i<n;i++){
			aux = 0;
			d=0;
			for (int k=nbits-1;k>=0;k--){
				fmpz_poly_get_coeff_mpz(tmp, poly_nums[i],k);
				d = fhe_decrypt(tmp, sk);
				aux = (aux*2)+d;
			} 
		printf("%d ", aux);
		
	}
	printf("\n");


	}
	
	////////////// Decryption /////////////////////


/*		for(i=0;i<n;i++){
		printf(" \n The %d'th number is \n", i);
		for(k=0;k<nbits;k++){
			fmpz_poly_get_coeff_mpz(tmp, poly_nums[i],k);
			printf(" %i ", fhe_decrypt(tmp, sk));
		}
	
	}
	
	printf("\n");	
	mpz_t tmp;	
	mpz_init(tmp);

	int d ;
	printf("\n After Bitonic sort \n");
	for (i=0;i<n;i++){
		aux = 0;
		d=0;
		for (int k=nbits-1;k>=0;k--){
			fmpz_poly_get_coeff_mpz(tmp, poly_nums[i],k);
			d = fhe_decrypt(tmp, sk);
			aux = (aux*2)+d;
		} 
		printf(" %d ", aux);
		
	}
*/	
	printf("\n");
	for(int k=0;k<n;k++)
		fmpz_poly_clear( poly_nums[k]);
	free(poly_nums);
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);	
	mpz_clear(c0);
	mpz_clear(tmp);
	free(list);			


}

void test_insertion_sort(){   // Generate n random numbers of nbits each and sort them
	int n = 32;// Nombre d'entiers à trier
	int nbits = 4; // Size in number of bits
	int *list = malloc(sizeof(int)*n); // List d'entiers à trier
	//printf("\n");

	mpz_t * max;
	mpz_t * min;
	max = malloc(sizeof(mpz_t) * nbits);
	min = malloc(sizeof(mpz_t) * nbits);
	for(int i=0;i<nbits;i++){
		mpz_init(max[i]);
		mpz_init(min[i]);
	}
	
	mpz_t tmp;	
	mpz_init(tmp);
	
	struct timeval start, end;
	long mtime, seconds, useconds;    
	

	fmpz_poly_t * poly_nums;
	poly_nums = malloc(sizeof(fmpz_poly_t) * n);
	
	int aux ; 	
	int a ;

	int i;
	int k;
	mpz_t c0;
	
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	fhe_keygen(pk, sk);
	
	mpz_init(c0);
	
	for(int obser = 0 ; obser < 3 ;obser++){	
		srand(time(NULL));
		int mod = (int)pow(2, nbits);
		//printf(" The array is \n");
	
		for(int t=0;t<n;t++){
			list[t]= rand()%mod;
		//printf("%d ", list[t]);
		}
	
	nbits=0;
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
			if(i+1>nbits) nbits = i+1;
		}

	//printf("Number of bits in this case is %d \n ", nbits);
int k;	


	/////////// Evaluation ////////////////////
	
	
	
		
	clock_t  START_eval = clock();
	gettimeofday(&start, NULL);

	
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

	double T_Elapsed4 = (double) (clock () - START_eval);
	//printf(" Evaluation took %f clock/sec \n ", T_Elapsed4);
					
	gettimeofday(&end, NULL);
	seconds  = end.tv_sec  - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
   	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
	

   	printf("%d \t %d \t %f \t %ld \n",n, nbits, T_Elapsed4, mtime);
  }
	
	//////////////// Decryption ///////////////////////////
/*	
	int d ;
	printf("\n After insertion sort \n");
	for (i=0;i<n;i++){
		aux = 0;
		d=0;
		for (k=nbits-1;k>=0;k--){
			fmpz_poly_get_coeff_mpz(tmp, poly_nums[i],k);
			d = fhe_decrypt(tmp, sk);
			aux = (aux*2)+d;
		} 
		printf(" %d ", aux);
		
	}
	
	printf("\n");
*/
	for(k=0;k<nbits;k++){
		mpz_clear(max[k]);
		mpz_clear(min[k]);
	}
	free(max);
	free(min);	
	for(k=0;k<n;k++)
		fmpz_poly_clear( poly_nums[k]);
	
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);	
	mpz_clear(c0);
	mpz_clear(tmp);
	free(list);
	
}

void test_sum_mpz_t(fmpz_poly_t sum_result, fmpz_poly_t poly_c1, fmpz_poly_t poly_c2, fhe_pk_t pk , fhe_sk_t sk ){	
	/////////// Evaluation ////////////////////
	int nbits,ninf, n1, n2;
	n1=fmpz_poly_degree(poly_c1)+1 ;
	n2=fmpz_poly_degree(poly_c2)+1 ;
	if(n1 > n2){
		nbits = n1; ninf = n2;
	} else {
		nbits = n2; ninf = n1;
	}
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

	mpz_t a_k, b_k, c_in, c_out, aux;
	mpz_init(a_k);	mpz_init(b_k);	mpz_init(c_in);	
	mpz_init(c_out); mpz_init(aux);
	fhe_encrypt(c_in, pk, 0);fhe_encrypt(aux, pk, 0);
	//fhe_encrypt(c_out, pk, 0);

	int k;
	for(k=0;k<ninf;k++){
		fmpz_poly_get_coeff_mpz(a_k, poly_c1,k);	
		fmpz_poly_get_coeff_mpz(b_k, poly_c2,k);
//		printf("for k = %d,  a_k is %i ", k, fhe_decrypt(a_k,sk));
//		printf(" b_k is %i  \n ", fhe_decrypt(b_k, sk));
		//printf("b's bit is %i \n", fhe_decrypt(b_k, sk));
		fhe_fulladd(tmp, c_out, a_k, b_k, c_in, pk);
	//	fmpz_poly_get_coeff_mpz(a_k, sum_result, k);
//		printf("sum at %d th loop is %i \n", k, fhe_decrypt(tmp, sk));
//		printf("Carry at %d th loop is %i\n",k,fhe_decrypt(c_out, sk));
		fmpz_poly_set_coeff_mpz(sum_result,k,tmp);
		mpz_set(c_in, c_out);					
	}
	if(n1 < n2){
		for(k=ninf;k<nbits;k++){
			fmpz_poly_get_coeff_mpz(b_k, poly_c2,k);

			fhe_fulladd(tmp, c_out, aux, b_k, c_in, pk);
		
			fmpz_poly_set_coeff_mpz(sum_result,k,tmp);
			mpz_set(c_in, c_out);					
		}
	} else {
		for(k=ninf;k<nbits;k++){
			fmpz_poly_get_coeff_mpz(a_k, poly_c1,k);	

			fhe_fulladd(tmp, c_out, a_k, aux, c_in, pk);

			fmpz_poly_set_coeff_mpz(sum_result,k,tmp);
			mpz_set(c_in, c_out);	
		}
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
	
	
	mpz_clear(c_in); mpz_clear(c_out); mpz_clear(tmp);
	mpz_clear(a_k); mpz_clear(b_k); mpz_clear(aux); 
}

void test_majority_bit(){
	
	unsigned aux1, aux2;
	int i, nbits ;

	mpz_t c0;
	mpz_t c1;
	fmpz_poly_t poly_c1;
	fmpz_poly_t poly_c2;

	mpz_t tmp;
	mpz_init(tmp);

	fmpz_poly_t tmp1_poly;
	fmpz_poly_init(tmp1_poly);
	fmpz_poly_t tmp2_poly;
	fmpz_poly_init(tmp2_poly);

	mpz_t tmp1 ;
	mpz_init(tmp1);
	mpz_t tmp2 ;
	mpz_init(tmp2);	

	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	fhe_keygen(pk, sk);
	
	mpz_init(c0);
	mpz_init(c1);
	fmpz_poly_init(poly_c1);
	fmpz_poly_init(poly_c2);
	struct timeval start, end;	
    	long mtime, seconds, useconds;    
	for(unsigned a = 45; a < 55; a++){
		printf("a: %d , ", a);
		gettimeofday(&start, NULL);   
		aux1= a * 50; 
		i=0;
		////// Encryption of the bit sequences ////////////////
		//printf("aux1[i] -----> %i\n", aux1 % 2);
		fhe_encrypt(c0, pk, a % 2);
		fmpz_poly_set_coeff_mpz( poly_c1 , i , c0 );
		i++;		
		aux1 = aux1 >> 1;
		do {		
			//printf("aux1[i] -----> %i\n", aux1 % 2);
			fhe_encrypt(c0, pk, aux1 % 2);
			fmpz_poly_set_coeff_mpz (poly_c1 , i , c0 );
			i++;
			aux1 = aux1 >> 1;	
		}while(aux1 != 0 );
		////////////////////////////////////////

		nbits=i;
		aux2=nbits/2;
	
		////// Encryption of the bit sequences ////////////////
		//printf("aux2[i] -----> %i\n", aux2 % 2);
		fhe_encrypt(c1, pk, nbits % 2);
		i=0;
		fmpz_poly_set_coeff_mpz(poly_c2, i, c1);		
		aux2 = aux2 >> 1;
		do {		
			//	printf("aux2[i]------>%i\n", aux2 % 2);
			fhe_encrypt(c1, pk, aux2 %2);
			i++;
			fmpz_poly_set_coeff_mpz (poly_c2, i, c1);
			aux2 = aux2 >> 1;	
		}while(aux2!=0 );
		////////////////////////////////////////	

		////////////////// Evaluation ///////////////////////	

		fhe_encrypt(tmp1, pk, 0);
		fmpz_poly_set_coeff_mpz(tmp1_poly,0,tmp1);
		//fmpz_poly_set_coeff_mpz(tmp2_poly,0,tmp1);
		//fhe_encrypt(tmp1,pk,1);
		//fmpz_poly_set_coeff_mpz(tmp2_poly,1,tmp1);
		//test_sum_mpz_t(tmp1_poly, tmp2_poly, tmp1_poly,pk,sk);
	
	
		for(long k = 0 ; k<nbits; k++){		
			fmpz_poly_get_coeff_mpz(tmp2, poly_c1, k);
			fmpz_poly_set_coeff_mpz(tmp2_poly,0,tmp2);

			//printf("\n tmp2_poly \n");
			//for (int j=fmpz_poly_degree(tmp2_poly);j>=0;j--){
			//	fmpz_poly_get_coeff_mpz(tmp1,tmp2_poly,j);
			//	printf("%i*", fhe_decrypt(tmp1,sk));
			//}
			//printf("\n tmp1_poly \n");
			//for (int j=fmpz_poly_degree(tmp1_poly);j>=0;j--){
			//	fmpz_poly_get_coeff_mpz(tmp1,tmp1_poly,j);
			//	printf("%i.", fhe_decrypt(tmp1,sk));
			//}		

			test_sum_mpz_t(tmp1_poly, tmp1_poly, tmp2_poly,pk, sk);

			//printf("\n sum_poly \n");
			//for (int j=fmpz_poly_degree(tmp1_poly);j>=0;j--){
			//fmpz_poly_get_coeff_mpz(tmp1,tmp1_poly,j);
			//printf("%i-", fhe_decrypt(tmp1,sk));
			//}			
		}
	
		//printf(" \n Printing the sum of n bits \n");

		/*  for(int k=fmpz_poly_degree(tmp1_poly);k>=0;k--){
		  fmpz_poly_get_coeff_mpz(tmp1,tmp1_poly,k);
		  printf(" %i ",fhe_decrypt(tmp1,sk));
		  }
		  printf(" \n");
		*/
		//// Calling isGreaterfunction //////////

		test_aIsGreater(tmp1,tmp1_poly, poly_c2, pk, fmpz_poly_degree(tmp1_poly)+1);

		//printf("The majority bit in %d is %i", a, fhe_decrypt(tmp1, sk) ) ;
		
		//printf("\n");
		gettimeofday(&end, NULL);
		seconds  = end.tv_sec  - start.tv_sec;
		useconds = end.tv_usec - start.tv_usec;
		mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
		printf("nbits: %d , temps: %ld \n",nbits,  mtime );
	}
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

void sum32_ui(fmpz_poly_t s, unsigned a , unsigned b, fhe_pk_t pk  , fhe_sk_t sk){
	unsigned aux1 = a, aux2 = b;
	int i=0;
	mpz_t c1, c2;
	mpz_init(c1); mpz_init(c2); 
	fmpz_poly_t poly_a, poly_b;
	fmpz_poly_init(poly_a); fmpz_poly_init(poly_b);
	//
	//printf("aux1[i] -----> %i\t", aux1 % 2);
	//printf("aux2[i] -----> %i\n", aux2 % 2);

	fhe_encrypt(c1, pk, aux1 % 2);
	fhe_encrypt(c2, pk, aux2 % 2);
	fmpz_poly_set_coeff_mpz( poly_a , i , c1 );
	fmpz_poly_set_coeff_mpz( poly_b , i , c2 );
	aux1 = aux1 >> 1;
	aux2 = aux2 >> 1;
	//
	while( aux1 != 0 || aux2 != 0){		
		//printf("aux1[i] -----> %i\t", aux1 % 2);
		//printf("aux2[i] -----> %i\n", aux2 % 2);
		fhe_encrypt(c1, pk, aux1 % 2);
		fhe_encrypt(c2, pk, aux2 % 2);
		i++;
		fmpz_poly_set_coeff_mpz (poly_a , i , c1 );
		fmpz_poly_set_coeff_mpz (poly_b , i , c2 );
		aux1 = aux1 >> 1;	
		aux2 = aux2 >> 1;	
	}
	//
	test_sum_mpz_t(s, poly_a, poly_b, pk , sk );
	printf("nb de bits %d ", i+1);
	//
	mpz_clear(c1); mpz_clear(c2); 
	fmpz_poly_clear(poly_a); fmpz_poly_clear(poly_b);
}
/*
void test_sum_unbounded(){
	mpz_t a ,b, res;
	mpz_init(a); mpz_init(b); mpz_init(res);
	///
	int i= mpz_set_str (a, "151315" , 0);
	assert(i==0);
	i= mpz_set_str (b, "105090601030102030405060708090" , 0);
	assert(i==0);
	///
	char * str_a = mpz_get_str(NULL, 10 , a); 
	char * str_b = mpz_get_str(NULL, 10 , b);

	printf("la chaine: %s de longueur %d \n", str_a, (int)strlen(str_a));
	printf("la chaine: %s de longueur %d \n", str_b,(int) strlen(str_b));

	char * s_ai;
	char * s_bi[40];
	strncpy(s_bi, str_b, 10);
	s_bi[10]= '\0';
	printf("la chaine: %s de longueur %d \n", s_bi, strlen(s_bi));

	//unsigned ai, bi;

	///
	mpz_clear(a); mpz_clear(b); mpz_clear(res);
}
*/

void mesure_sum_integers(){
	//KeyGen
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	fhe_keygen(pk, sk);

	struct timeval start, end;	
	long mtime, seconds, useconds; 
   
	unsigned a = 0, b= 1;
	for(int i= 0; i < 30; i++){
		a= a*2; b = (2 * b )/2;
		for(int k= 1; k < 100; k++){
			a= a + k ; b = b + 3 * k ;

			gettimeofday(&start, NULL);   

			test_sum_integers(a , b, pk, sk);

			gettimeofday(&end, NULL);
			seconds  = end.tv_sec  - start.tv_sec;
			useconds = end.tv_usec - start.tv_usec;
			mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
			printf("time: %ld ms \n", mtime );
		}
	}
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);

}

void test_sum_integers(unsigned a , unsigned b, fhe_pk_t pk, fhe_sk_t sk){	

	unsigned res;
	fmpz_poly_t s;
	//
	mpz_t c;

	fmpz_poly_init(s);
	mpz_init(c);
	 // Integers to be added
	printf("a: %u , b: %u ", a,b);
	//
	sum32_ui(s, a, b, pk, sk);	
	//decrypt
	res=0;
	for(int i=fmpz_poly_degree(s); i >= 0; i--){
		fmpz_poly_get_coeff_mpz(c,s,i);
		res = res * 2 + fhe_decrypt(c, sk);
	}
	printf("sum: %u ", res); assert(res == a + b);
	//
	fmpz_poly_clear(s);
	mpz_clear(c);

}

void mesure_test_matrix_prod(){
	/*
	int m, n, p, q, c, d, k ;
  	
	int first[10][10], second[10][10];
 	mpz_t first_enc[10][10];
	mpz_t second_enc[10][10];
	mpz_t multiply_enc[10][10];
	*/
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	fhe_keygen(pk, sk);
	struct timeval start, end;	
    	long mtime, seconds, useconds;    
	//gettimeofday(&start, NULL);  
	for(int i= 4; i <7; i++){
		gettimeofday(&start, NULL);
		test_matrix_prod( 4*i,4* i, 4*i,4* i,pk, sk);
		gettimeofday(&end, NULL);
		seconds  = end.tv_sec  - start.tv_sec;
		useconds = end.tv_usec - start.tv_usec;
		mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
		printf("size : %d , time : %ld ms \n",4*i, mtime );
	}
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);
}
void test_matrix_prod(int m, int n, int p, int q,fhe_pk_t pk,fhe_sk_t sk){
	
	int c, d, k ;
  	
	int first[30][30], second[30][30];
 	mpz_t first_enc[30][30];
	mpz_t second_enc[30][30];
	mpz_t multiply_enc[30][30];
	

	mpz_t sum;
	mpz_init(sum);
	mpz_t tmp;
	mpz_init(tmp);
	
	/*

	printf("Enter the number of rows and columns of first matrix\n");
  	scanf("%d%d", &m, &n);
	printf("Enter the elements of first matrix\n");
 
 	for ( c = 0 ; c < m ; c++ ){
   		for ( d = 0 ; d < n ; d++ ){
     			scanf("%d", &first[c][d]);
			mpz_init(first_enc[c][d]);
			fhe_encrypt(first_enc[c][d], pk, first[c][d]);
			
		}
 	}

	
 	 printf("Enter the number of rows and columns of second matrix\n");
 	 scanf("%d%d", &p, &q);
	*/
  	if ( n != p )
 	   	printf("Matrices with entered orders can't be multiplied with each other.\n");
	else{
			
   		printf("Enter the elements of second matrix\n");
 		fhe_encrypt(sum, pk, 0);	 
		for ( c = 0 ; c < p ; c++ ){
      			for ( d = 0 ; d < q ; d++ ){
       				scanf("%d", &second[c][d]);
 				mpz_init(second_enc[c][d]);
				fhe_encrypt(second_enc[c][d], pk, second[c][d]);
			}
		}


   		 for ( c = 0 ; c < m ; c++ ){
      			for ( d = 0 ; d < q ; d++ ){
      			 	 for ( k = 0 ; k < p ; k++ ){
					fhe_mul(tmp, first_enc[c][k], second_enc[k][d], pk);
					fhe_add(sum, sum, tmp, pk);
         				
       				 }
 	
        			mpz_init(multiply_enc[c][d]);
				mpz_set(multiply_enc[c][d],sum);
        			fhe_encrypt(sum, pk,0);
     			 }
    		}
 
 	   	printf("Product of entered matrices:-\n");
 
   		 for ( c = 0 ; c < m ; c++ ) {
      			for ( d = 0 ; d < q ; d++ )
        			printf("%d\t", fhe_decrypt(multiply_enc[c][d], sk));
 
      			printf("\n");
   		 }
  	}	



}


void test_xor_bits(){
	unsigned m, aux;
	mpz_t c0, c1;
	fmpz_poly_t poly_c;
	mpz_init(c0);
	mpz_init(c1);
	fmpz_poly_init(poly_c);
	fhe_pk_t pk;
	fhe_sk_t sk;
	fhe_pk_init(pk);
	fhe_sk_init(sk);
	fhe_keygen(pk, sk);
	int j, nbits;	

	struct timeval start, end;	
	long mtime, seconds, useconds;    
	//clock_t  START_eval;
	//double T_Elapsed4;

	m= 2147483648	;
	for(int i=0; i< 250; i++){
		m= m + 2;
		j=0;
		//printf("--------->%i\n", m % 2);
		fhe_encrypt(c0, pk, m % 2);
		fmpz_poly_set_coeff_mpz ( poly_c , 0 , c0 );
		aux = m;	
		aux = aux >> 1;
		do {
			//	printf("--------->%i\n", aux % 2);
			fhe_encrypt(c0, pk, aux % 2);
			aux = aux >> 1;
			j++;
			fmpz_poly_set_coeff_mpz ( poly_c , j , c0 );
			//fhe_add(c0, c0, c1, pk);
		}while(aux != 0);
		nbits= j+1;
		//START_eval = clock();
		gettimeofday(&start, NULL);    
		mpz_set(c1, c0);
		while(j>0){
			j--;
			fmpz_poly_get_coeff_mpz ( c0 , poly_c , j );
			fhe_add(c1, c0, c1, pk);
		}
		gettimeofday(&end, NULL);
		seconds  = end.tv_sec  - start.tv_sec;
		useconds = end.tv_usec - start.tv_usec;
		mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
		//printf("Elapsed time in KeyGen : %ld ms \n", mtime );

		aux = fhe_decrypt(c1, sk);
		printf("xor de %u: %d bits, runtime: %ld ms\n", m, nbits, mtime);
	}
	mpz_clear(c0);
	mpz_clear(c1);
	fmpz_poly_clear(poly_c);
	fhe_pk_clear(pk);
	fhe_sk_clear(sk);
}

void test_encryt_decrypt(){
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
		fhe_encrypt(c0, pk, 0);
		//fhe_encrypt(c1, pk, 1);
		//printf("\nadd-chain: ");
		//for (int j = 0; j < RUNS*RUNS; j++) {
		//	fhe_add(c0, c0, c0, pk);
		//	fhe_recrypt(c0, pk);
		//	m = fhe_decrypt(c0, sk);
		//	printf("%i", m);
		//	fflush(stdout);
		//}
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
