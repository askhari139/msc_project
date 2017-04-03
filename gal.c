#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)



float ran2(long *idum);// function declaration for random number generation
float prod_rate(float a, float b, float c, float d);//function declaration for production rate
float rate(float a, int b, int c);//for other rates to make sure that the rates arezero in case of non-positive conc.
void main()
{

	FILE *output, *parameters;// define the files
	
	int i,conc[7],j,k, count, z;// define integer variables: conc[] for concentrations, i,j and k loop variables
	float prod[8], degrade[7], inter[6], rates[17], prob[17], sum_prob[17], time, x, gal, fact, total_rate, hill[3];//define float variables: parameters, rates and probabilities
	x=1.0;
	gal=x/100;	
	int loop=0;
	for (loop =0; loop <10; loop++)
	{
	parameters = fopen ("vent_para.txt", "r"); //open the parmeters file with read permissions
	for (i=0; i<7; i++)// read concentrations: gal1*, gal3*, gal4, gal80, complexes, gal1 and gal3 respectively
		{
			fscanf(parameters, "%d", &conc[i]);
		}
	for (i=0; i<8; i++)
		{fscanf (parameters, "%f", &prod[i]);}
	for (i=0; i<7; i++) //results change when this is set as 4
		{fscanf (parameters, "%f", &degrade[i]);}
	for(i=0; i<6; i++)
		{fscanf (parameters, "%f", &inter[i]);}
	fscanf(parameters, "%f", &fact);
	for (i=0; i<3; i++)
	{
		fscanf(parameters, "%f", &hill[i]);
	}
	fclose(parameters);
	long id = -7432483288.0+loop+x; 
	time = 0;
	count = 0;// time recording variable
	char s[sizeof"gal_10000.txt"];
	sprintf(s, "gal_%d.txt", loop);
	output = fopen(s, "w");
	while (time < 20000)
	{
		for (i = 0; i < 17; i++)
		{
			sum_prob[i] = 0;
		}
		float choice = ran2(&id);
		float t = ran2(&id);
		//--------------------------defining the rates-------------------------
		//---------------------------------------------------------------------
		//production rates
		rates[0] = prod_rate(prod[0], prod[1], conc[2], hill[0])+ gal*fact;// gal1
		rates[1] = prod_rate(prod[2], prod[3], conc[2], hill[1]) + gal;// gal3
		rates[2] = prod[6];//gal 4
		rates[3] = prod[7]+ prod_rate(prod[4], prod[5], conc[2], hill[2]);// gal80
		//degradation rates
		for (i=4; i<11; i++)
			{rates [i] = rate(degrade[i-4],conc[i-4],1);}

		//protein reaction rates
		rates[11] = rate(inter[0],conc[0],conc[3]); //gal1 + gal80
		rates[12] = rate(inter[1],conc[4],1);
		rates[13] = rate(inter[2],conc[1],conc[3]);// gal3 + gal80
		rates[14] = rate(inter[3],conc[5],1);
		rates[15] = rate(inter[4],conc[2],conc[3]);// gal4 +gal80
		rates[16] = rate(inter[5],conc[6],1);
		//--------------------------------------------------------------------------

		//total rate

		total_rate = 0;
		for (i=0; i<17; i++)
			{total_rate += rates[i]; }
		
		//--------------------------------------------------------------------------

		//probabilities
		for (i = 0; i < 17; ++i)
		{
			prob[i] = rates[i]/total_rate;
		}
		
		for (i=0; i<17; i++)
		{
			for (j =0; j<=i; j++)
			{
				sum_prob[i] += prob[j];
			}
			
		}
			if (choice<= sum_prob[0])
				{conc[0] +=1;}
				
			else if (choice<=sum_prob[1])
				{conc[1] += 1;}
				
			else if (choice<=sum_prob[2])
				{conc[2] += 1;}
				
			else if (choice<=sum_prob[3])
				{conc[3] += 1;}
		//-----------------------------------------------------------------------------

		//degradation------------------------------------------------------------------		
			else if (choice<=sum_prob[4])
				{conc[0] -= 1;}
				
			else if (choice<=sum_prob[5])
				{conc[1] -= 1;}
				
			else if (choice<=sum_prob[6])
				{conc[2] -= 1;}
				
			else if (choice<=sum_prob[7])
				{conc[3] -=1;}
			else if (choice<=sum_prob[8])
				{conc[4] -= 1;}
				
			else if (choice<=sum_prob[9])
				{conc[5] -= 1;}
				
			else if (choice<=sum_prob[10])
				{conc[6] -= 1;}
				
		//--------------------------------------------------------------------------------

		//gal80 reaction------------------------------------------------------------------
			else if (choice<=sum_prob[11])
				{conc[0] -= 1;
				conc[3] -= 1;
				conc[4] += 1;}
				
			else if (choice<=sum_prob[12])
				{conc[0] += 1;
				conc[3] += 1;
				conc[4] -= 1;}
				
			else if (choice<=sum_prob[13])
				{conc[1] -= 1;
				conc[3] -= 1;
				conc[5] += 1;}
				
			else if (choice<=sum_prob[14])
				{conc[1] += 1;
				conc[3] += 1;
				conc[5] -= 1;}
				
			else if (choice<=sum_prob[15])
				{conc[2] -= 1;
				conc[3] -= 1;
				conc[6] += 1;}
				
			else
				{conc[2] += 1;
				conc[3] += 1;
				conc[6] -= 1;}

		time += -1.0*log(t)/total_rate;
		z = (int)(time);
		//if (z >= count)
		//{ 
			fprintf(output, "%f %d %d %d %d %d %d %d\n", time, conc[0], conc[1], conc[2], conc[3], conc[4], conc[5], conc[6]);
			count += 1;
		//}
		
	}
	fclose(output);
	}

}


float prod_rate(float a, float b, float c, float d)
{
	float rate, power;
	power = pow(c,b);
	rate = a * (power/(pow(d,b) + power));
	if (c >0)
	{return rate;}
	else
	{return 0;}
}

float rate(float a, int b, int c)
{
	float r = a*b*c;
	if (b<=0 || c <=0)
		{return 0;}
	else
		{return r;}
}


float ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

 