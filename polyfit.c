/*Functtion :多项式拟合polyfit**********************************************/
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef unsigned char    boolean;
typedef int polifie_mth;enum
    {
    LAT_IS_VARIABLE,
    LON_IS_VARIABLE,
    TIME_IS_VARIABLE,
    };

typedef int derection_mode;enum
    {
    DERECTION_1,
    DERECTION_2,
    DERECTION_3,
    DERECTION_4,
    DERECTION_5,
    DERECTION_6,
    DERECTION_7,
    DERECTION_8,
    };

typedef struct 
    {
    double averange_rate;
    double delta_lat;
    double delta_lon;
    }track_direction_intf;

#define SUM_DATA     10000
#define DATE_LENGTH  100
#define STR_LINE_LENGTH  1024
#define AVG_NUM     3


#define AREA_STR "Area: "
#define LENGTH_STR "Length: "
#define CIRFCE_STR "Circumference: "
#define ANGLE_STR "Angle: "

#define TRUE 1
#define FALSE -1
int main(int argc,char *argv[])
{
    int i,j,m,n=7,poly_n=2;
    int txt_line_num;
    double x[7]={1,2,3,4,6,7,8},y[7]={2,3,6,7,5,3,2};
    double a[3];
    
    
    double lat[ SUM_DATA ]={0},lon[SUM_DATA]={0},ele[SUM_DATA]={0};
    double time[SUM_DATA] = {0};
    track_direction_intf area_track_inf;     
    
    double lat_max = 0.0f;
    double lat_min = 0.0f;
    double lon_max = 0.0f;
    double lon_min = 0.0f;
    double k_lat[3];
    double k_lon[3];
    double k_averange;

    char date_string[SUM_DATA][DATE_LENGTH] = {'\0'};
    char string_line[STR_LINE_LENGTH];
    char *data_start;
    char *data_offset;
    char *data_end;
    char data_string[20];
    int max_min_is_init = FALSE;


    void polyfit(int n,double *x,double *y,int poly_n,double a[]);
    track_direction_intf get_track_rate_info(  int n,double *lon,double *lat );
    void insertsort_positive(double array[],int len);
    void insertsort_invert(double array[],int len);
    void insertsort(double array[],int len,boolean positive);

    //const char file_path[] = {"C:\\Users\\zhaobruce\\Desktop\\fix_length\\20181029\\eTrex302\\Area_2018-10-26_182104.txt"};
    char file_path[100];
    //const char file_path_new[] ={"C:\\Users\\zhaobruce\\Desktop\\fix_length\\20181015\\Area_2018-10-15_100902_new.txt"};
    char file_path_new[100];// ={"C:\\Users\\zhaobruce\\Desktop\\fix_length\\20181015\\Area_2018-10-15_100902_new.txt"};
    system("cls");
    for(i=0;i<argc;i++)
    {
        printf("i=%d\n",i);
        printf(argv[i]);
        printf("\n");
    }
    if( argc < 2 )
       {
        return -1;
       }
    //strcpy(argv[1],"C:\\Users\\zhaobruce\\Desktop\\fix_length\\20181015\\Area_2018-10-15_100902.txt"); 
    //strcpy(file_path,"C:\\Users\\zhaobruce\\Desktop\\fix_length\\20181015\\Area_2018-10-15_100902.txt");
    //strcpy(data_start,"C:\\Users\\zhaobruce\\Desktop\\fix_length\\20181015\\Area_2018-10-15_100902.txt");
    strcpy(file_path,argv[1]);
    data_offset = strstr(file_path,".txt");
    strncpy(file_path_new, file_path, data_offset - file_path );
    strcat(file_path_new,"_new.txt");
    FILE *fp_read = fopen( file_path,"r");
    if( NULL == fp_read )
        {
            printf("File Path Error");
            return 0;
        }

    FILE *fp_write = fopen( file_path_new,"w");
    if(fp_write==NULL)  
        {              
        printf("File Path Error");
        return 0;  
        }

    i = 0;
    fgets( string_line,STR_LINE_LENGTH,fp_read );
    
    if(strstr(string_line, AREA_STR ) != NULL)
    {
        fprintf(fp_write,"%s",string_line);
        fgets( string_line,STR_LINE_LENGTH,fp_read );
    }

    if(strstr(string_line, LENGTH_STR ) != NULL)
    {
        fprintf(fp_write,"%s",string_line);
        fgets( string_line,STR_LINE_LENGTH,fp_read );
    }

    if(strstr(string_line, CIRFCE_STR ) != NULL)
    {
        fprintf(fp_write,"%s",string_line);
        fgets( string_line,STR_LINE_LENGTH,fp_read );
    }

    if(strstr(string_line, ANGLE_STR ) != NULL)
    {
        fprintf(fp_write,"%s",string_line);
        fgets( string_line,STR_LINE_LENGTH,fp_read );
    }
    //data_offset = -strlen(string_line);
   // max_min_is_init = fseek(fp_read,data_offset,2); 
    if( NULL != string_line )
        {
        i = 0;
        sscanf( string_line,"%s\t%lf\t%lf\t%lf",date_string[i],&lat[i],&lon[i],&ele[i] );
        #if(0)
        const char file_path[] = {"C:\\Users\\zhaobruce\\Desktop\\fix_length\\20181015\\Area_2018-10-15_100902.txt"}; = string_line;
        data_offset = strstr(string_line,"\t");
        strncpy(date_string[i],string_line,data_offset - data_start );

        data_start = data_offset + 1;
        data_offset = strstr(data_start,"\t");
        strncpy(data_string,data_start,data_offset - data_start );
        data_string[data_offset - data_start + 1] = '\0';
        lat[i] = atof( data_string );
        lat_max = lat[i];
        lat_min = lat[i];

        data_start = data_offset + 1;
        data_offset = strstr(data_start,"\t");
        strncpy(data_string,data_start,data_offset - data_start );
        data_string[data_offset - data_start + 1] = '\0';
        lon[i] = atof( data_string );
        lon_max = lon[i];
        lon_min = lon[i];

        data_start = data_offset + 1;
        data_offset = strstr(data_start,"\n");
        strncpy(data_string,data_start,data_offset - data_start );
        data_string[data_offset - data_start ] = '\0';
        ele[i] = atof( data_string );
        #endif
        i++;
        }

    while( fscanf( fp_read,"%s\t%lf\t%lf\t%lf",date_string[i],&lat[i],&lon[i],&ele[i] ) != EOF )
        {
        if(lat[i] > lat_max )
            {
            lat_max = lat[i];
            }

        if(lat[i] < lat_min )
            {
            lat_min = lat[i];
            }
        
        if(lon[i] > lon_max )
            {
            lon_max = lon[i];
            }
            
        if(lon[i] < lon_min )
            {
            lon_min = lon[i];
            }
        //printf("%s\t%.6f\t%.6f\t%.6f\n",date_string[i],lat[i],lon[i],ele[i]);
        time[i] = i + 1;
        i++;
        }
    fclose(fp_read);
    txt_line_num = i;
    //polyfit(n,x,y,poly_n,a);
    //polyfit(i,lat,lon,poly_n,a);
   area_track_inf = get_track_rate_info( txt_line_num, lon, lat);
    printf("k_averange=%f\n",area_track_inf.averange_rate);
   if( ( area_track_inf.averange_rate >= -1 ) && ( area_track_inf.averange_rate <= 1 ) )
        {
        polyfit(txt_line_num,lon,lat,poly_n,k_lon);
        for (i=0;i<poly_n+1;i++)/*这里是升序排列，Matlab是降序排列*/
            {
            printf("k_lon[%d]=%.10f\n",i,k_lon[i]);
            }

        if( area_track_inf.delta_lon > 0 )
            {
            insertsort_positive( lon, txt_line_num );
            }
        else
            {
            insertsort_invert( lon, txt_line_num );
            }
        
        for( i = 0; i < txt_line_num; i++)
            {
            lat[i] = k_lon[0] + k_lon[1] * lon[i] + k_lon[2] * lon[i] * lon[i];
            fprintf(fp_write,"%s\t%.6f\t%.6f\t%.6f\n",date_string[i],lat[i],lon[i],ele[i]);
            //printf("%s\t%.6f\t%.6f\t%.6f\n",date_string[i],lat[i],lon[i],ele[i]);            
            }
        }
    else
        {
        polyfit(txt_line_num,lat,lon,poly_n,k_lat);
        for (i=0;i<poly_n+1;i++)/*这里是升序排列，Matlab是降序排列*/
            {
            printf("k_lat[%d]=%.10f\n",i,k_lat[i]);
            }

        if( area_track_inf.delta_lat > 0 )
            {
            insertsort( lat, txt_line_num ,1 );
            }
        else
            {
            insertsort( lat, txt_line_num, 0 );
            }
        
        for( i = 0; i < txt_line_num; i++)
            {
            lon[i] = k_lat[0] + k_lat[1] * lat[i] + k_lat[2] * lat[i] *  lat[i];
            //lon[i] = k_lon[0] + k_lon[1] * time[i] + k_lon[2] * time[i] *  time[i];
            fprintf(fp_write,"%s\t%.6f\t%.6f\t%.6f\n",date_string[i],lat[i],lon[i],ele[i]);
            //printf("%s\t%.6f\t%.6f\t%.6f\n",date_string[i],lat[i],lon[i],ele[i]);            
            }
       }
/*
   {
        default:
            polyfit(txt_line_num,time,lat,poly_n,k_lat);
            polyfit(txt_line_num,time,lon,poly_n,k_lon);
            for (i=0;i<poly_n+1;i++)/*这里是升序排列，Matlab是降序排列*/
 /*               printf("k_lat[%d]=%g\n",i,k_lat[i]);
            for (i=0;i<poly_n+1;i++)/*这里是升序排列，Matlab是降序排列*/
  /*              printf("k_lon[%d]=%g\n",i,k_lon[i]);

            for( i = 0; i < txt_line_num; i++)
                {
                lat[i] = k_lat[0] + k_lat[1] * time[i] + k_lat[2] * time[i] *  time[i];
                lon[i] = k_lon[0] + k_lon[1] * time[i] + k_lon[2] * time[i] *  time[i];
                if( ( lat[i] <= lat_max ) && ( lat[i] >= lat_min ) && ( lon[i] <= lon_max ) && ( lon[i] >= lon_min ) )
                        {
                        fprintf(fp_write,"%s\t%.6f\t%.6f\t%.6f\n",date_string[i],lat[i],lon[i],ele[i]);
                        //printf("%s\t%.6f\t%.6f\t%.6f\n",date_string[i],lat[i],lon[i],ele[i]);            
                        }
                    else
                        {
                        printf("%s\t%.6f\t%.6f\t%.6f\n",date_string[i],lat[i],lon[i],ele[i]);  
                        }
                }
            printf("lat_max=%.6f\nlat_min=%.6f\nlon_max=%.6f\nlon_min=%.6f\n",lat_max,lat_min,lon_max,lon_min);
            break;
   }
*/
    fclose(fp_write);
    printf("Sucessful!\n"); 
    return 0;

    getch();
}

/*==================polyfit(n,x,y,poly_n,a)===================*/
/*=======拟合y=a0+a1*x+a2*x^2+……+apoly_n*x^poly_n========*/
/*=====n是数据个数 xy是数据值 poly_n是多项式的项数======*/
/*===返回a0,a1,a2,……a[poly_n]，系数比项数多一（常数项）=====*/
void polyfit(int n,double x[],double y[],int poly_n,double a[])
{
    int i,j;
    double *tempx,*tempy,*sumxx,*sumxy,*ata;
    void gauss_solve(int n,double A[],double x[],double b[]);
    tempx=(double*)calloc(n,sizeof(double));
    sumxx=(double*)calloc(poly_n*2+1,sizeof(double));
    tempy=(double*)calloc(n,sizeof(double));
    sumxy=(double*)calloc(poly_n+1,sizeof(double));
    ata=(double*)calloc((poly_n+1)*(poly_n+1),sizeof(double));
    for (i=0;i<n;i++)
        {
        tempx[i]=1;
        tempy[i]=y[i];
        }
    for (i=0;i<2*poly_n+1;i++)
        {
        for (sumxx[i]=0,j=0;j<n;j++)
            {
            sumxx[i]+=tempx[j];
            tempx[j]*=x[j];
            }
        }
    for (i=0;i<poly_n+1;i++)
        {
        for (sumxy[i]=0,j=0;j<n;j++)
            {
            sumxy[i]+=tempy[j];
            tempy[j]*=x[j];
            }
        }

    for (i=0;i<poly_n+1;i++)
        {
        for (j=0;j<poly_n+1;j++)
            {
            ata[i*(poly_n+1)+j]=sumxx[i+j];
            }
        }
    gauss_solve(poly_n+1,ata,a,sumxy);
    free(tempx);
    free(sumxx);
    free(tempy);
    free(sumxy);
    free(ata);
}

void gauss_solve(int n,double A[],double x[],double b[])
{
    int i,j,k,r;
    double max;
    for(k=0;k<n-1;k++)
        {
        max=fabs(A[k*n+k]); /*find maxmum*/
        r=k;
        for(i=k+1;i<n-1;i++)
            {
            if(max<fabs(A[i*n+i]))
                {
                max=fabs(A[i*n+i]);
                r=i;
                }
            }
        if(r!=k)
        for (i=0;i<n;i++)   /*change array:A[k]&A[r] */
            {
            max=A[k*n+i];
            A[k*n+i]=A[r*n+i];
            A[r*n+i]=max;
            }
        max=b[k];/*change array:b[k]&b[r] */
        b[k]=b[r];
        b[r]=max;
        for(i=k+1;i<n;i++)
            {
            for(j=k+1;j<n;j++)
            A[i*n+j]-=A[i*n+k]*A[k*n+j]/A[k*n+k];
            b[i]-=A[i*n+k]*b[k]/A[k*n+k];
            }
        }
    for(i=n-1;i>=0;x[i]/=A[i*n+i],i--)
    for(j=i+1,x[i]=b[i];j<n;j++)
    x[i]-=A[i*n+j]*x[j];
}

track_direction_intf get_track_rate_info( int n,double *lon,double *lat )
{
    int i,j;
    int average_num = AVG_NUM;
    double *tempx,*tempy,*sumxx,*sumxy,*ata;
    double averange_lat_start;
    double averange_lat_end;
    double averange_lon_start;
    double averange_lon_end;
    double delta_lat;
    double delta_lon;
    double averange_k;
    track_direction_intf track_infor;
    //tempx=(double*)calloc(n,sizeof(double));
    //tempy=(double*)calloc(n,sizeof(double));
    for( i=0; i < AVG_NUM; i++ )
        {
        averange_lat_start += lat[i];
        averange_lon_start += lon[i];
        averange_lat_end   += lat[n - i - 1];
        averange_lon_end   += lon[n - i - 1];
        }

    averange_lat_start = averange_lat_start / 3;
    averange_lon_start = averange_lon_start / 3;
    averange_lat_end   = averange_lat_end / 3;
    averange_lon_end   = averange_lon_end / 3;

    delta_lat = averange_lat_end - averange_lat_start;
    delta_lon = averange_lon_end - averange_lon_start;
    averange_k = delta_lat / delta_lon;
    track_infor.averange_rate = averange_k;
    track_infor.delta_lat     = delta_lat;
    track_infor.delta_lon     = delta_lon;
    return  track_infor;
}

void insertsort_positive(double array[],int len)
{
	int i, j;
    double temp;
	for (i = 1; i < len; i++) {
		temp = array[i];
		j = i - 1;
		while (j >= 0 && array[j] > temp) {
			array[j + 1] = array[j];
			j--;
		}
		if (j != i - 1)
			array[j + 1] = temp;
	}
}

void insertsort_invert(double array[],int len)
{
	int i, j;
    double temp;
	for (i = 1; i < len; i++) {
		temp = array[i];
		j = i - 1;
		while (j >= 0 && array[j] < temp) {
			array[j + 1] = array[j];
			j--;
		}
		if (j != i - 1)
			array[j + 1] = temp;
	}
}

void insertsort(double array[],int len,boolean positive)
{
int i, j;
double temp;
for (i = 1; i < len; i++)
    {
    temp = array[i];
    j = i - 1;
    if( positive )
        {
        while (j >= 0 && array[j] > temp)
            {
            array[j + 1] = array[j];
            j--;
            }
        }
    else
        {
        while (j >= 0 && array[j] < temp)
            {
            array[j + 1] = array[j];
            j--;
            }
        }
    if (j != i - 1)
    array[j + 1] = temp;
    }
}