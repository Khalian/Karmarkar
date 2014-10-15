#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
 
uint32_t         stampstart();
uint32_t         stampstop(uint32_t start);
int reset(float x[50],int tracker,int n,float arrlower[50]); 
void check(float constants[50][50],float x[50],float b[50],int eqno,int n, uint32_t start);
void mul(float a[50][50],float b[50][50],int n,int m, int p,int q);

void display(float x[50],int n);
int reset(float x[50],int tracker,int n,float arrlower[50]);

int main()
{
	uint32_t         start, stop;

int i,j,tracker,n,eqno;                       //total nos declaration still to be done
float arrlower[50],arrupper[50],x[50],tempval,totalno=1;
float constants[50][50],b[50];

printf("\n Enter the dimension of the problem:");
scanf("%d",&n);

printf("\n Enter the lower bounds:");
for(i=0;i<n;i++)
scanf("%f",&arrlower[i]);

printf("\n Enter the upper bounds:");
for(i=0;i<n;i++)
scanf("%f",&arrupper[i]);

printf("\n Enter the no of equations:");
scanf("%d",&eqno);

printf("\n Enter the elements of the constraint matrix:");
for(i=0;i<eqno;i++)
for(j=0;j<n;j++)
scanf("%f",&constants[i][j]);

printf("\n Enter the elements of the b matrix:");
for(i=0;i<eqno;i++)
scanf("%f",&b[i]);

printf("\n");
	start = stampstart();


for(i=0;i<n;i++)
x[i]=arrlower[i];        //lowest possible value of x coord/vector initialised   

tracker=n-1;             // the initial value of the indice that tracks the position of the coord to increment

for(i=0;i<n;i++)
{
tempval=arrupper[i]-arrlower[i]+1;
totalno=totalno*tempval;
}

printf("\n Worst case no of iterations is:%f",totalno);

for(i=0;i<totalno;i++)         //this is the principal loop
{

if(x[tracker]>=arrlower[tracker]&&x[tracker]<=arrupper[tracker])    //base case to check at trackr position
{
check(constants,x,b,eqno,n,start);
x[tracker]++;
}

else
{  //start of else
i--;
here:
tracker--;

if(x[tracker]>=arrlower[tracker]&&x[tracker]<=(arrupper[tracker]-1))
{
x[tracker]++;
tracker=reset(x,tracker,n,arrlower);
}

else
goto here;                                 
}  //end of else

}                              //end of principal loop 

printf("\n There does not exist any integer solutions");

printf("\n");
	stop = stampstop(start);
}                              //end of main program

int reset(float x[50],int tracker,int n,float arrlower[50])
{
int i;
for(i=(tracker+1);i<n;i++)     //bounds need reviewing
x[i]=arrlower[i];

return(n-1);

}
void check(float constants[50][50],float x[50],float b[50],int eqno,int n,uint32_t start)
{
    int i,j,flag=1;                 //the initial assumtion that its true
    float x1[50][50],a[50][50];
    uint32_t stop;

    for(i=0;i<n;i++)
    {
        x1[i][0]=x[i];
    }

    for(i=0;i<eqno;i++)
        for(j=0;j<n;j++)
            a[i][j]=constants[i][j];

    mul(a,x1,eqno,n,n,1);
    
    for(i=0;i<eqno;i++)
    {
        if(a[i][0]!=b[i])
        flag=0;                        //reassigment if at all it is false
    }
    if(flag==1)
    {
        printf("\n An integer solution exists,it occurs at:");
        for(i=0;i<n;i++)
            printf("%f,",x[i]);
        stop = stampstop(start);
        exit(0);
    }
}

void mul(float a[50][50],float b[50][50],int n,int m, int p,int q)    //mul is replacing a with a*b
{
    float c[50][50];
    int i,j,k;

    if(p!=m)
    {
        printf("\n Matrix multiplication not possible");
        exit(0);
    }

    for(i=0;i<n;i++)
        for(j=0;j<q;j++)
            c[i][j]=0;

    for(i=0;i<n;i++)
        for(j=0;j<q;j++)
            for(k=0;k<m;k++)
                c[i][j]=c[i][j]+a[i][k]*b[k][j];

    for(i=0;i<n;i++)
        for(j=0;j<q;j++)
            a[i][j]=c[i][j];
}

uint32_t
stampstart()
{
	struct timeval  tv;
	struct timezone tz;
	struct tm      *tm;
	uint32_t         start;
 
	gettimeofday(&tv, &tz);
	tm = localtime(&tv.tv_sec);
 
	printf("TIMESTAMP-START\t  %d:%02d:%02d:%d (~%d ms)\n", tm->tm_hour,
	       tm->tm_min, tm->tm_sec, tv.tv_usec,
	       tm->tm_hour * 3600 * 1000 + tm->tm_min * 60 * 1000 +
	       tm->tm_sec * 1000 + tv.tv_usec / 1000);
 
	start = tm->tm_hour * 3600 * 1000 + tm->tm_min * 60 * 1000 +
		tm->tm_sec * 1000 + tv.tv_usec / 1000;
 
	return (start);
 
}
 
uint32_t
stampstop(uint32_t start)
{
 
	struct timeval  tv;
	struct timezone tz;
	struct tm      *tm;
	uint32_t         stop;
 
	gettimeofday(&tv, &tz);
	tm = localtime(&tv.tv_sec);
 
	stop = tm->tm_hour * 3600 * 1000 + tm->tm_min * 60 * 1000 +
		tm->tm_sec * 1000 + tv.tv_usec / 1000;
 
	printf("TIMESTAMP-END\t  %d:%02d:%02d:%d (~%d ms) \n", tm->tm_hour,
	       tm->tm_min, tm->tm_sec, tv.tv_usec,
	       tm->tm_hour * 3600 * 1000 + tm->tm_min * 60 * 1000 +
	       tm->tm_sec * 1000 + tv.tv_usec / 1000);
 
	printf("ELAPSED\t  %d ms\n", stop - start);
 
	return (stop);
 
}

