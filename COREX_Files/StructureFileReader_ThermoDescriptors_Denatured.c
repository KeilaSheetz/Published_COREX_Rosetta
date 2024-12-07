/**** Structure File Reader Application ****/
/**** Calculates Residue-Specific Thermo Descriptors ****/
/**** Changed int2bin ****/

/*** STW 9/2002, RSS ***/

#define         Str_Len             100
#define         Default_FName       "PDB Filename"
#define         Default_Answer      "Yes"
#define         Max_Atoms           5000
#define         Num_AA	            30
#define         Max_Res             600
#define         Max_WindowSize      30
#define         Max_Nunits          100
#define         Def_window_size     5
#define         Def_Min_Res_Win     4
#define         Def_Nfrac           10

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <stdlib.h>
typedef char bstr[Str_Len];


/* Global symbols */

float *Atom_Area_Native;                               /* main */
bstr inFName;
char atomname[Max_Atoms][6];
char resname[Max_Res][8];
int First_Res_Num;
int Last_Res_Num;
int window_size;
int Min_Res_Win;
double Prob_unfolded25[Max_Res];
double Partition_Function25;
int partition;
float Fraction_Folded;
float Sconf;
float delASA_ap;
float delASA_pol;
FILE *outfile;
int pdb_resnum[Max_Res];
char pdb_molecule[Max_Res][6];
double deltaH_pol[Max_Res];
double deltaH_ap[Max_Res];
double deltaSconf[Max_Res];
double deltaG_solv[Max_Res];
double H_ratio[Max_Res];
double S_ratio[Max_Res];
double k_f[Max_Res];
double lnk_f[Max_Res];
char code[Max_Res][6];
float env_odds[Max_Res];

char *aa3[Num_AA];                                     /* Ecsact_Par */
float aCp;
float adeltaH;
float ASAexapol[Num_AA];
float ASAexOH[Num_AA];
float ASAexpol[Num_AA];
float ASAsc[Num_AA];
float bCp;
float bdeltaH;
float dSbb[Num_AA];
float dSbbGly;
float dSbuex[Num_AA];
float dSexu[Num_AA];
float dSbb_length_corr;
float minusAla;
float minusBeta;
float minusOther;
float mw[Num_AA];
float OHCp;
float TsApolar;
float TsPolar;
float W_Sconf;
float Current_Temp;
float exposed_criteria;
float ASA_exposed_amide;
int NAATotal;
char *thermo_env[13];
float log_odds[Num_AA][13];

FILE *corefile;                                        /* LogFileOpen */

int Npartition;                                      /* Gen_Partition */
int First_Residue[Max_WindowSize][Max_Nunits];
int Last_Residue[Max_WindowSize][Max_Nunits];
int Nunits[Max_WindowSize];

char stateflag[Max_Nunits];                            /* int2bin */

double dG25;                                           /* calc_component_energies */
double stat_weight;
float dH_ap;
float dH_pol;
float dSconf;
float dG_solv;

double dH_pol_nf[Max_Res];                              /* Residue_Probabilities */
double dH_ap_nf[Max_Res];
double dSconf_nf[Max_Res];
double dG_solv_nf[Max_Res];
double dH_pol_f[Max_Res];
double dH_ap_f[Max_Res];
double dSconf_f[Max_Res];
double dG_solv_f[Max_Res];


/* Subroutines */

extern void Ecsact_Par (void);
extern void LogFileOpen (void);
extern void Gen_Partition (void);
extern void int2bin (long k);
extern void calc_component_energies (void);
extern void Residue_Probabilities (void);
extern void Save_Probabilities (void);

void ZotNewLine();

/************************/

int main (int argc,char**argv)

{
int npointsmax;
bstr instring;
int num_atoms;
int Total_Res;
int oldresnum;
int OTnum;
int atomnum;
int resnum;
int i;
char label_nam[7];
float A_Area=0;
char name[7];
int Begin_Atom_Res[Max_Res];
int Last_Atom_Res[Max_Res];
int j;
bstr out_fname;
bstr inFName2;
FILE *infile;
float temp;
float best; /* JOW 5-16-22 */
char Answer[10];
long state_count;
long state_number;
int temp_num;
char temp_char[6];
double Part_Funct_nf[Max_Res];
double Part_Funct_f[Max_Res];
int k;

npointsmax=Max_Atoms;
Atom_Area_Native=(float*)malloc(npointsmax*sizeof(float));

Ecsact_Par();

printf("\nEnter the original PDB file name [%s] > ",Default_FName);
    fgets(inFName,Str_Len,stdin);
    ZotNewLine(inFName);
    if (!strlen(inFName))
        strcpy(inFName,Default_FName);

LogFileOpen();

fgets(instring,Str_Len,corefile);
sscanf(instring,"%d",&num_atoms);
fgets(instring,Str_Len,corefile);
Total_Res=0;
oldresnum=-100;
OTnum=0;
atomnum=0;
resnum=0;
for(i=0;i<num_atoms;i++)
   {
    fgets(instring,Str_Len,corefile);
    sscanf(instring,"%3c%d%d%4c%f%6d%3s",name,&resnum,&atomnum,label_nam,&A_Area,&temp_num,temp_char);
    name[4]='\0';
    label_nam[4] = '\0';
    temp_char[4]='\0';
    pdb_resnum[resnum]=temp_num;
    strncpy(pdb_molecule[resnum],temp_char,4);
    Atom_Area_Native[atomnum]=A_Area;
    if (resnum!=oldresnum)
       {
        strncpy(resname[resnum],name,3);
        Begin_Atom_Res[resnum]=atomnum;
        Last_Atom_Res[resnum-1]=atomnum-1;
       }
    strncpy(atomname[atomnum],label_nam,4);
    if (strncmp(atomname[atomnum]," OXT",4)==0||strncmp(atomname[atomnum]," OT",3)==0) OTnum++;
    if (strncmp(atomname[atomnum]," CA",3)==0) Total_Res++;
    oldresnum=resnum;
    if(atomnum==0) First_Res_Num=resnum;
   }
Last_Atom_Res[resnum]=atomnum;
Last_Res_Num=First_Res_Num+Total_Res-1;
printf("Total Number of residues = %d. Last atom %s %d\n",Total_Res,atomname[atomnum],atomnum);
printf("First Residue %s %d %s, Last Residue %s %d %s\n",resname[First_Res_Num],pdb_resnum[First_Res_Num],pdb_molecule[First_Res_Num],resname[Last_Res_Num],pdb_resnum[Last_Res_Num],pdb_molecule[Last_Res_Num]);
fclose(corefile);

printf("\nEnter Entropy Scaling Factor [%5.3f] > ",W_Sconf); 
    fgets(instring,Str_Len,stdin);
    ZotNewLine(instring);
    if (!strlen(instring))
        W_Sconf=1.5;
    else
       W_Sconf=1.5;

printf("\nEnter Temperature [%2.0f deg C] > ",Current_Temp-273.15); 
    fgets(instring,Str_Len,stdin);
    ZotNewLine(instring);
    if (!strlen(instring))
        Current_Temp=298.15;
    else
       {
        sscanf(instring,"%f",&Current_Temp);
        Current_Temp = Current_Temp+273.15;
       }

printf("\nEnter Size of Folding Units [%d] > ",Def_window_size); 
    fgets(instring,Str_Len,stdin);
    ZotNewLine(instring);
    if (!strlen(instring))
        window_size=Def_window_size;
    else
        sscanf(instring,"%d",&window_size);

printf("\nEnter Minimum Size of Folding Unit [%d] > ",Def_Min_Res_Win); 
    fgets(instring,Str_Len,stdin);
    ZotNewLine(instring);
    if (!strlen(instring))
        Min_Res_Win=Def_Min_Res_Win;
    else
       sscanf(instring,"%d",&Min_Res_Win);

Gen_Partition();

for(i=First_Res_Num;i<=Last_Res_Num;i++)
   {
    Prob_unfolded25[i]=0.0;

    /* I added these initializations in 7/20/09 JOW */
    dH_pol_nf[i]=0.0;
    dH_ap_nf[i]=0.0;
    dSconf_nf[i]=0.0;
    dG_solv_nf[i]=0.0;
    dH_pol_f[i]=0.0;
    dH_ap_f[i]=0.0;
    dSconf_f[i]=0.0;
    dG_solv_f[i]=0.0;
    deltaH_pol[i]=0.0;
    deltaH_ap[i]=0.0;
    deltaSconf[i]=0.0;
    deltaG_solv[i]=0.0;
    H_ratio[i]=0.0;
    S_ratio[i]=0.0;
    k_f[i]=0.0;
    lnk_f[i]=0.0;
    /* end of JOW initializations */
   }
Partition_Function25=1.0;

printf("\nMonte Carlo Sampled Structure Files (Y/N) [%s] > ",Default_Answer);
    fgets(instring,Str_Len,stdin);
    ZotNewLine(instring);
    if (!strlen(instring))
        strcpy(Answer,Default_Answer);
    else
        strcpy(Answer,instring);
if ((strncmp(Answer,"Y",1)==0)||(strncmp(Answer,"y",1)==0))
   {
    sprintf(inFName,"%s.MC",inFName);
   }
sprintf(inFName2,"%s.%d.%d",inFName,window_size,Min_Res_Win);
if ((infile=fopen(inFName2,"r"))==NULL) 
   {
    fprintf(stderr,"\nfopen() error #%d on structure file %s\n",errno,inFName2);
    exit(1);
   }

printf("\nRunning... ");
printf("\n");
while(fscanf(infile,"%d",&partition))
     {
      if(feof(infile))
	break;
      fscanf(infile,"%ld",&state_count);
      fscanf(infile,"%ld",&state_number);
      fscanf(infile,"%f",&Fraction_Folded);
      fscanf(infile,"%f",&Sconf);
      fscanf(infile,"%f",&delASA_ap);
      fscanf(infile,"%f",&delASA_pol);
      int2bin(state_number);
      printf("\r");
      /*printf("partition %5d    state number %6ld     stateflag  %s",partition,state_number,stateflag);*/
      printf("partition %5d    state number %6ld",partition,state_number);
      /*printf("partition %5d",partition);*/
      fflush(stdout);
      for(i=1;i<=Nunits[partition];i++)
         {
          if (stateflag[i-1]=='0')
             {
              for(j=First_Residue[partition][i];j<=Last_Residue[partition][i];j++)
                 {
                  fscanf(infile,"%f",&temp);
                 }
             }
         }
      while (getc(infile)!='\n');
      calc_component_energies();
      Residue_Probabilities();
     }
fclose(infile);

for(i=First_Res_Num;i<=Last_Res_Num;i++)
   {
    Part_Funct_nf[i]=Prob_unfolded25[i];
    Part_Funct_f[i]=Partition_Function25-Part_Funct_nf[i];
    deltaH_pol[i]=dH_pol_f[i]/Part_Funct_f[i]-dH_pol_nf[i]/Part_Funct_nf[i];
    deltaH_ap[i]=dH_ap_f[i]/Part_Funct_f[i]-dH_ap_nf[i]/Part_Funct_nf[i];
    deltaSconf[i]=dSconf_f[i]/Part_Funct_f[i]-dSconf_nf[i]/Part_Funct_nf[i];
    deltaG_solv[i]=dG_solv_f[i]/Part_Funct_f[i]-dG_solv_nf[i]/Part_Funct_nf[i];
    Prob_unfolded25[i]=Prob_unfolded25[i]/Partition_Function25;
   }

for(i=First_Res_Num;i<=Last_Res_Num;i++)
   {
    H_ratio[i]=deltaH_pol[i]/deltaH_ap[i];
    S_ratio[i]=deltaSconf[i]/deltaG_solv[i];
    k_f[i]=(1.0-Prob_unfolded25[i])/Prob_unfolded25[i];
    lnk_f[i]=log(k_f[i]);
   }

/* Here is where TE are assigned from Descriptors, if S==0.500 or S==1.500, JOW 5-16-22 */
if(W_Sconf==1.5){
 for(i=First_Res_Num;i<=Last_Res_Num;i++)
   {
    best=1000000.0; code[i][0]='0'; code[i][1]='0'; code[i][2]='0'; code[i][3]='\0';
    temp=abs(-1.9872*Current_Temp*lnk_f[i]-10696.0)+abs(deltaH_ap[i]-(-1538.0))+abs(deltaH_pol[i]-1094.0)+abs(Current_Temp*deltaSconf[i]-(-9939.0));
    if(temp<best) { code[i][0]='1'; code[i][1]='1'; code[i][2]='1'; code[i][3]='\0'; best=temp; }
    temp=abs(-1.9872*Current_Temp*lnk_f[i]-10178.0)+abs(deltaH_ap[i]-(-1522.0))+abs(deltaH_pol[i]-2426.0)+abs(Current_Temp*deltaSconf[i]-(-8363.0));
    if(temp<best) { code[i][0]='2'; code[i][1]='2'; code[i][2]='2'; code[i][3]='\0'; best=temp; }
    temp=abs(-1.9872*Current_Temp*lnk_f[i]-9859.0)+abs(deltaH_ap[i]-(-1028.0))+abs(deltaH_pol[i]-(-1075.0))+abs(Current_Temp*deltaSconf[i]-(-10902.0));
    if(temp<best) { code[i][0]='3'; code[i][1]='3'; code[i][2]='3'; code[i][3]='\0'; best=temp; }
    temp=abs(-1.9872*Current_Temp*lnk_f[i]-9100.0)+abs(deltaH_ap[i]-662.0)+abs(deltaH_pol[i]-(-400.0))+abs(Current_Temp*deltaSconf[i]-(-10327.0));
    if(temp<best) { code[i][0]='4'; code[i][1]='4'; code[i][2]='4'; code[i][3]='\0'; best=temp; }
    temp=abs(-1.9872*Current_Temp*lnk_f[i]-9082.0)+abs(deltaH_ap[i]-(-996.0))+abs(deltaH_pol[i]-44.0)+abs(Current_Temp*deltaSconf[i]-(-9312.0));
    if(temp<best) { code[i][0]='5'; code[i][1]='5'; code[i][2]='5'; code[i][3]='\0'; best=temp; }
    temp=abs(-1.9872*Current_Temp*lnk_f[i]-8427.0)+abs(deltaH_ap[i]-(-1115.0))+abs(deltaH_pol[i]-1085.0)+abs(Current_Temp*deltaSconf[i]-(-7704.0));
    if(temp<best) { code[i][0]='6'; code[i][1]='6'; code[i][2]='6'; code[i][3]='\0'; best=temp; }
    temp=abs(-1.9872*Current_Temp*lnk_f[i]-7611.0)+abs(deltaH_ap[i]-(-303.0))+abs(deltaH_pol[i]-(-2168.0))+abs(Current_Temp*deltaSconf[i]-(-9790.0));
    if(temp<best) { code[i][0]='7'; code[i][1]='7'; code[i][2]='7'; code[i][3]='\0'; best=temp; }
    temp=abs(-1.9872*Current_Temp*lnk_f[i]-7500.0)+abs(deltaH_ap[i]-(-1083.0))+abs(deltaH_pol[i]-(-915.0))+abs(Current_Temp*deltaSconf[i]-(-8323.0));
    if(temp<best) { code[i][0]='8'; code[i][1]='8'; code[i][2]='8'; code[i][3]='\0'; best=temp; }
    /* Following is the old code */
    /*code[i][0]='M';
    if (lnk_f[i]<7.95) code[i][0]='L';
    if (lnk_f[i]>=13.4) code[i][0]='H';
    temp=-1.024*deltaH_ap[i]-2553.0;
    code[i][1]='L';
    if (deltaH_pol[i]>temp) code[i][1]='H';
    temp=0.125*deltaG_solv[i]-3053.0;
    code[i][2]='L';
    if ((Current_Temp*deltaSconf[i])>temp) code[i][2]='H';
    code[i][3]='\0';*/
   }
}
else { for(i=First_Res_Num;i<=Last_Res_Num;i++){ code[i][0]='0'; code[i][1]='0'; code[i][2]='0'; code[i][3]='\0'; } }
/* Here is where scoring table is encoded, from the array at tail-end, JOW 5-16-22 */
for(i=First_Res_Num;i<=Last_Res_Num;i++) env_odds[i]=0.0;

if(W_Sconf==1.5){
 for(i=First_Res_Num;i<=Last_Res_Num;i++)
   {
    j=0;
    k=0;
    for(j=1;j<NAATotal;j++)
       {
        if (strncmp(resname[i],aa3[j],3)==0)
           {
            for(k=1;k<=8;k++) /* k used to run to 12, JOW 5-16-22 */
               {
                if (strcmp(code[i],thermo_env[k])==0)
                   {
                    env_odds[i]=log_odds[j][k];
		    code[i][1]='\0';
                    break;
                   }
               }
           }
       }
   }
}

sprintf(out_fname,"%s.W%d.T%2.0fS%2.3f.ThermoDescriptD",inFName,window_size,Current_Temp-273.15,W_Sconf);
printf("\nOut file = %s\n",out_fname);
outfile=fopen(out_fname,"w");
Save_Probabilities();
exit(1);

}

/******************** LogFileOpen for Partition Files  **************/

void LogFileOpen()
{
char FileNull,Q1;
bstr temp;

FileNull='N';
Q1='N';
while(FileNull=='N')
     {
      if (Q1=='Y')
         {
          printf("\nEnter name of original PDB file: ");
          fgets(inFName,Str_Len,stdin);
          ZotNewLine(inFName);
          sprintf(temp,"%s.info",inFName);
         }
      else
         {
          sprintf(temp,"%s.info",inFName);
         }
      corefile=fopen(temp,"r");
      if (corefile==NULL)
         {
          printf("\nA file named %s.info does not exist in the current directory.",inFName);
          Q1='Y';
         }
      else
         {
          FileNull='Y';
         }
     }

}

/************************** Gen_Partition ******************/

void Gen_Partition (void)
{
int i;
int partition_index;
int end;

partition_index=1;
Npartition=window_size;
Last_Residue[1][0]=First_Res_Num-1;
for(partition_index=1;partition_index<=Npartition;partition_index++)
   {
    end=1;
    First_Residue[partition_index][1]=First_Res_Num;
    Last_Residue[partition_index][1]=First_Res_Num+partition_index-2;
    i=1;
    if (partition_index==1) i=0;
    while(end)
         {
          i++;
          First_Residue[partition_index][i]=Last_Residue[partition_index][i-1]+1;
          Last_Residue[partition_index][i]=First_Residue[partition_index][i]+window_size-1;
          if (Last_Residue[partition_index][i]>=Last_Res_Num)
             {
              Last_Residue[partition_index][i]=Last_Res_Num;
              end=0;
              Nunits[partition_index]=i;
             }
         }
   }

              /***  This adjusts the partition to allow for a minimum number of residues per window ****/

for(partition_index=1;partition_index<=Npartition;partition_index++)
   {
    if (Last_Residue[partition_index][1]-First_Residue[partition_index][1]<Min_Res_Win-1)
       {
        Last_Residue[partition_index][1]=Last_Residue[partition_index][2];
        for(i=2;i<Nunits[partition_index];i++)
           {
            First_Residue[partition_index][i]=First_Residue[partition_index][i+1];
            Last_Residue[partition_index][i]=Last_Residue[partition_index][i+1];
           }
        Nunits[partition_index]=Nunits[partition_index]-1;
       }
    if (Last_Residue[partition_index][Nunits[partition_index]]-First_Residue[partition_index][Nunits[partition_index]]<Min_Res_Win-1)
       {
        Last_Residue[partition_index][Nunits[partition_index]-1]=Last_Residue[partition_index][Nunits[partition_index]];
        Nunits[partition_index]=Nunits[partition_index]-1;
       }
   }

}

/************************************ int2bin *****************************************/

/* Converts integer state number to binary code */

void int2bin (long k)			
{
char s[Max_Nunits];
long i;

strcpy(stateflag,"0000");
for(i=0;i<Nunits[partition];i++) s[i]='0';
s[Nunits[partition]]='\0';
for(i=0;i<(sizeof(k))*8;i++)
   {
    if ((k>>i)&1) s[i]='1';
   }
sprintf(stateflag,"%s",s);

}

/******************************** calc_component_energies *********************************/

void calc_component_energies(void)
{
float SconfN,dSsolv_ap,dSsolv_pol;

SconfN=Sconf*W_Sconf;

dH_ap=delASA_ap*(adeltaH+aCp*((Current_Temp-273.15)-60));
dH_pol=delASA_pol*(bdeltaH+bCp*((Current_Temp-273.15)-60));
dSconf=SconfN;
dSsolv_ap=delASA_ap*aCp*log((Current_Temp/TsApolar));
dSsolv_pol=delASA_pol*bCp*log((Current_Temp/TsPolar));
dG_solv=dH_ap+dH_pol-Current_Temp*(dSsolv_ap+dSsolv_pol);
dG25=dG_solv-Current_Temp*SconfN;
stat_weight=exp(-dG25/(1.9872*Current_Temp));
Partition_Function25=Partition_Function25+stat_weight;

}

/******************************* Residue_Probabilities ***************************/

void Residue_Probabilities()
{
int i,j;

for(i=1;i<=Nunits[partition];i++)
   {
    if (stateflag[i-1]=='1')
       {
        for(j=First_Residue[partition][i];j<=Last_Residue[partition][i];j++)
           {
            Prob_unfolded25[j]=Prob_unfolded25[j]+stat_weight;
            dH_pol_nf[j]=dH_pol_nf[j]+dH_pol*stat_weight;
            dH_ap_nf[j]=dH_ap_nf[j]+dH_ap*stat_weight;
            dSconf_nf[j]=dSconf_nf[j]+dSconf*stat_weight;
            dG_solv_nf[j]=dG_solv_nf[j]+dG_solv*stat_weight;
           }
       }
    if (stateflag[i-1]=='0')
       {
        for(j=First_Residue[partition][i];j<=Last_Residue[partition][i];j++)
           {
            dH_pol_f[j]=dH_pol_f[j]+dH_pol*stat_weight;
            dH_ap_f[j]=dH_ap_f[j]+dH_ap*stat_weight;
            dSconf_f[j]=dSconf_f[j]+dSconf*stat_weight;
            dG_solv_f[j]=dG_solv_f[j]+dG_solv*stat_weight;
           }
       }
   }

}

/********************************** Save_Probabilities **************************/

void Save_Probabilities(void)
{
int i;

fprintf(outfile,"Res  Index  dHpol  dHap  dSconf  dGsolv  H_ratio  S_ratio  k_f  ln(k_f)  DTE  Odds\n");
for(i=First_Res_Num;i<=Last_Res_Num;i++)
   {
    fprintf(outfile,"%s %d %e %e %e %e %e %e %e %e %1s %e\n",resname[i],pdb_resnum[i],deltaH_pol[i],deltaH_ap[i],deltaSconf[i],deltaG_solv[i],H_ratio[i],S_ratio[i],k_f[i],lnk_f[i],code[i],env_odds[i]);
   }

fclose(outfile);

}

/*******************ZotNewLine *******************/

void ZotNewLine(char *theStr)
{

theStr[strcspn(theStr, "\n")] = '\0';

}

/*******************************************/

void Ecsact_Par (void)
{
int i;

/* AMINOACID CONSTANTS --------------------------------------------------- */
/* mw of water 18, must be subtracted per peptide bond */
/* dS are after Aquino */
/* Areas are Irene optimized */

dSbbGly=6.5;           /* original 4.8 */
minusAla=2.4;          /* original 2.4 */
minusBeta=4.32;         /* original 3.5 */
minusOther=3.1;         /* original 2.8 */
for (i=1;i<=Num_AA;i++)
	{
	ASAexOH[i]=0.0;
	}

aa3[1]="ALA";
mw[1]=89.09;
ASAexapol[1]=70.0;
ASAexpol[1]=36.13;
ASAsc[1]=52.12;
dSbuex[1]=0.0;
dSexu[1]=0.0;
dSbb[1]=dSbbGly-minusAla;

aa3[2]="ARG";
mw[2]=174.2;
ASAexapol[2]=87.13;
ASAexpol[2]=126.06;
ASAsc[2]=169.19;
dSbuex[2]=7.11;
dSexu[2]=-0.84;
dSbb[2]=dSbbGly-minusOther;

aa3[3]="ASN";
mw[3]=132.12;
ASAexapol[3]=38.14;
ASAexpol[3]=104.0;
ASAsc[3]=113.78;
dSbuex[3]=3.29;
dSexu[3]=2.24;
dSbb[3]=dSbbGly-minusOther;

aa3[4]="ASP";
mw[4]=133.1;
ASAexapol[4]=42.06;
ASAexpol[4]=95.0;
ASAsc[4]=102.43;
dSbuex[4]=2.0;
dSexu[4]=2.16;
dSbb[4]=dSbbGly-minusOther;

aa3[5]="ASX";
mw[5]=132.61;
ASAexapol[5]=42.06;
ASAexpol[5]=95.0;
ASAsc[5]=107.5;
dSbuex[5]=2.64;
dSexu[5]=2.2;
dSbb[5]=dSbbGly-minusOther;

aa3[6]="CYS";
mw[6]=121.15;
ASAexapol[6]=30.32;
ASAexpol[6]=75.05;
ASAsc[6]=91.93;
dSbuex[6]=3.55;
dSexu[6]=0.61;
dSbb[6]=dSbbGly-minusOther;

aa3[7]="GLN";
mw[7]=146.15;
ASAexapol[7]=65.0;
ASAexpol[7]=121.55;
ASAsc[7]=128.73;
dSbuex[7]=5.02;
dSexu[7]=2.12;
dSbb[7]=dSbbGly-minusOther;

aa3[8]="GLU";
mw[8]=147.13;
ASAexapol[8]=71.13;
ASAexpol[8]=94.33;
ASAsc[8]=117.54;
dSbuex[8]=3.53;
dSexu[8]=2.27;
dSbb[8]=dSbbGly-minusOther;

aa3[9]="GLX";
mw[9]=146.64;
ASAexapol[9]=71.13;
ASAexpol[9]=94.33;
ASAsc[9]=115.47;
dSbuex[9]=4.25;
dSexu[9]=2.2;
dSbb[9]=dSbbGly-minusOther;

aa3[10]="GLY";
mw[10]=75.07;
ASAexapol[10]=26.17;
ASAexpol[10]=43.12;
ASAsc[10]=0.2;
dSbuex[10]=0.0;
dSexu[10]=0.0;
dSbb[10]=dSbbGly;

aa3[11]="HIS";
mw[11]=155.16;
ASAexapol[11]=90.01;
ASAexpol[11]=68.0;
ASAsc[11]=144.29;
dSbuex[11]=3.44;
dSexu[11]=0.79;
dSbb[11]=dSbbGly-minusOther;

aa3[12]="ILE";
mw[12]=131.17;
ASAexapol[12]=110.68;
ASAexpol[12]=10.877;
ASAsc[12]=103.81;
dSbuex[12]=1.74;
dSexu[12]=0.67;
dSbb[12]=dSbbGly-minusBeta;

aa3[13]="LEU";
mw[13]=131.17;
ASAexapol[13]=122.34;
ASAexpol[13]=27.5;
ASAsc[13]=134.53;
dSbuex[13]=1.63;
dSexu[13]=0.25;
dSbb[13]=dSbbGly-minusOther;

aa3[14]="LYS";
mw[14]=146.19;
ASAexapol[14]=101.28;
ASAexpol[14]=79.0;
ASAsc[14]=156.94;
dSbuex[14]=5.86;
dSexu[14]=1.02;
dSbb[14]=dSbbGly-minusOther;

aa3[15]="MET";
mw[15]=149.21;
ASAexapol[15]=104.55;
ASAexpol[15]=64.0;
ASAsc[15]=158.01;
dSbuex[15]=4.55;
dSexu[15]=0.58;
dSbb[15]=dSbbGly-minusOther;

aa3[16]="PHE";
mw[16]=165.19;
ASAexapol[16]=186.81;
ASAexpol[16]=36.5;
ASAsc[16]=176.73;
dSbuex[16]=1.40;
dSexu[16]=1.51;
dSbb[16]=dSbbGly-minusOther;

aa3[17]="PRO";
mw[17]=115.13;
ASAexapol[17]=100.77;
ASAexpol[17]=15.61;
ASAsc[17]=67.3;
dSbuex[17]=0.0;
dSexu[17]=0.0;
dSbb[17]=dSbbGly-minusOther;

aa3[18]="SER";
mw[18]=105.09;
ASAexapol[18]=55.12;
ASAexpol[18]=81.94;
ASAsc[18]=68.58;
dSbuex[18]=3.68;
dSexu[18]=0.55;
dSbb[18]=dSbbGly-minusOther;
ASAexOH[18]=39.0;

aa3[19]="THR";
mw[19]=119.12;
ASAexapol[19]=79.46;
ASAexpol[19]=41.12;
ASAsc[19]=105.29;
dSbuex[19]=3.31;
dSexu[19]=0.48;
dSbb[19]=dSbbGly-minusOther;
ASAexOH[19]=35.0;

aa3[20]="TRP";
mw[20]=204.23;
ASAexapol[20]=184.52;
ASAexpol[20]=52.28;
ASAsc[20]=222.74;
dSbuex[20]=2.74;
dSexu[20]=1.15;
dSbb[20]=dSbbGly-minusOther;

aa3[21]="TYR";
mw[21]=181.19;
ASAexapol[21]=175.76;
ASAexpol[21]=71.1;
ASAsc[21]=188.46;
dSbuex[21]=2.78;
dSexu[21]=1.74;
dSbb[21]=dSbbGly-minusOther;

aa3[22]="VAL";
mw[22]=117.15;
ASAexapol[22]=88.69;
ASAexpol[22]=17.82;
ASAsc[22]=105.62;
dSbuex[22]=0.12;
dSexu[22]=1.29;
dSbb[22]=dSbbGly-minusBeta;

aa3[23]="IVA";
mw[23]=85.12;
ASAexapol[23]=167.85;
ASAexpol[23]=30.73;
ASAsc[23]=132.2;
dSbuex[23]=1.63;
dSexu[23]=0.25;
dSbb[23]=dSbbGly-minusBeta;

aa3[24]="STA";
mw[24]=157.15;
ASAexapol[24]=164.81;
ASAexpol[24]=62.63;
ASAsc[24]=117.82;
dSbuex[24]=1.63;
dSexu[24]=0.25;
dSbb[24]=dSbbGly-1.7;

aa3[25]="ACE";
mw[25]=43.3;
ASAexapol[25]=85.5;
ASAexpol[25]=26.1;
ASAsc[25]=78.7;
dSbuex[25]=0.0;
dSexu[25]=0.0;
dSbb[25]=0.0;

NAATotal=25;

/*thermo_env[1]="LHH";
thermo_env[2]="LHL";
thermo_env[3]="LLH";
thermo_env[4]="LLL";
thermo_env[5]="MHH";
thermo_env[6]="MHL";
thermo_env[7]="MLH";
thermo_env[8]="MLL";
thermo_env[9]="HHH";
thermo_env[10]="HHL";
thermo_env[11]="HLH";
thermo_env[12]="HLL";*/

thermo_env[1]="111";
thermo_env[2]="222";
thermo_env[3]="333";
thermo_env[4]="444";
thermo_env[5]="555";
thermo_env[6]="666";
thermo_env[7]="777";
thermo_env[8]="888";


/* ALA */
log_odds[1][1]=-0.46;
log_odds[1][2]=-0.08;
log_odds[1][3]=-0.55;
log_odds[1][4]=-0.30;
log_odds[1][5]=0.26;
log_odds[1][6]=0.19;
log_odds[1][7]=0.14;
log_odds[1][8]=0.46;
/*log_odds[1][1]=0.447798;
log_odds[1][2]=-0.27829;
log_odds[1][3]=0.421961;
log_odds[1][4]=-0.14661;
log_odds[1][5]=0.153051;
log_odds[1][6]=-0.46381;
log_odds[1][7]=0.249315;
log_odds[1][8]=-0.24133;
log_odds[1][9]=0.040228;
log_odds[1][10]=-0.58219;
log_odds[1][11]=0.2937;
log_odds[1][12]=-0.55706;*/
/* ARG */
log_odds[2][1]=-1.86;
log_odds[2][2]=-0.03;
log_odds[2][3]=-1.83;
log_odds[2][4]=-0.47;
log_odds[2][5]=-0.53;
log_odds[2][6]=0.94;
log_odds[2][7]=0.35;
log_odds[2][8]=0.80;
/*log_odds[2][1]=-0.43267;
log_odds[2][2]=-0.41432;
log_odds[2][3]=0.308753;
log_odds[2][4]=-0.47702;
log_odds[2][5]=-0.27789;
log_odds[2][6]=-1.74974;
log_odds[2][7]=-0.11669;
log_odds[2][8]=0.091232;
log_odds[2][9]=-0.09203;
log_odds[2][10]=-0.21946;
log_odds[2][11]=0.834623;
log_odds[2][12]=0.546607;*/
/* ASN */
log_odds[3][1]=0.36;
log_odds[3][2]=-0.60;
log_odds[3][3]=0.85;
log_odds[3][4]=-0.00;
log_odds[3][5]=-0.06;
log_odds[3][6]=-1.20;
log_odds[3][7]=-0.29;
log_odds[3][8]=-1.29;
/*log_odds[3][1]=-0.41075;
log_odds[3][2]=-0.41075;
log_odds[3][3]=-0.14194;
log_odds[3][4]=0.062687;
log_odds[3][5]=-0.54053;
log_odds[3][6]=-0.46215;
log_odds[3][7]=0.545618;
log_odds[3][8]=0.417881;
log_odds[3][9]=-0.8304;
log_odds[3][10]=-0.3729;
log_odds[3][11]=-0.02585;
log_odds[3][12]=0.680225;*/
/* ASP */
log_odds[4][1]=0.18;
log_odds[4][2]=-0.79;
log_odds[4][3]=0.73;
log_odds[4][4]=-0.11;
log_odds[4][5]=0.04;
log_odds[4][6]=-1.14;
log_odds[4][7]=0.16;
log_odds[4][8]=-0.55;
/*log_odds[4][1]=-0.44228;
log_odds[4][2]=-0.10808;
log_odds[4][3]=0.27852;
log_odds[4][4]=0.260001;
log_odds[4][5]=-0.40775;
log_odds[4][6]=-0.12174;
log_odds[4][7]=0.297621;
log_odds[4][8]=0.634039;
log_odds[4][9]=-0.54347;
log_odds[4][10]=-0.64559;
log_odds[4][11]=-0.17071;
log_odds[4][12]=-0.0102;*/
/* ASX */
log_odds[5][1]=0.0;
log_odds[5][2]=0.0;
log_odds[5][3]=0.0;
log_odds[5][4]=0.0;
log_odds[5][5]=0.0;
log_odds[5][6]=0.0;
log_odds[5][7]=0.0;
log_odds[5][8]=0.0;
log_odds[5][9]=0.0;
log_odds[5][10]=0.0;
log_odds[5][11]=0.0;
log_odds[5][12]=0.0;
/* CYS */
log_odds[6][1]=0.38;
log_odds[6][2]=0.75;
log_odds[6][3]=-0.58;
log_odds[6][4]=-0.56;
log_odds[6][5]=-0.06;
log_odds[6][6]=0.19;
log_odds[6][7]=-1.07;
log_odds[6][8]=-0.66;
/*log_odds[6][1]=-0.82883;
log_odds[6][2]=-0.62816;
log_odds[6][3]=-0.24157;
log_odds[6][4]=0.433061;
log_odds[6][5]=-0.68098;
log_odds[6][6]=-1.27044;
log_odds[6][7]=-0.13383;
log_odds[6][8]=0.570532;
log_odds[6][9]=-0.54704;
log_odds[6][10]=-0.3592;
log_odds[6][11]=0.143851;
log_odds[6][12]=0.920546;*/
/* GLN */
log_odds[7][1]=-0.38;
log_odds[7][2]=-1.11;
log_odds[7][3]=0.34;
log_odds[7][4]=0.33;
log_odds[7][5]=-0.04;
log_odds[7][6]=-0.53;
log_odds[7][7]=0.55;
log_odds[7][8]=0.11;
/*log_odds[7][1]=-1.35675;
log_odds[7][2]=-0.93294;
log_odds[7][3]=-0.20987;
log_odds[7][4]=-0.3025;
log_odds[7][5]=-0.98575;
log_odds[7][6]=-0.39656;
log_odds[7][7]=0.25454;
log_odds[7][8]=0.181198;
log_odds[7][9]=-0.20519;
log_odds[7][10]=0.33455;
log_odds[7][11]=0.586289;
log_odds[7][12]=0.76992;*/
/* GLU */
log_odds[8][1]=0.10;
log_odds[8][2]=-0.51;
log_odds[8][3]=0.48;
log_odds[8][4]=0.29;
log_odds[8][5]=0.06;
log_odds[8][6]=-0.85;
log_odds[8][7]=0.12;
log_odds[8][8]=-0.54;
/*log_odds[8][1]=-0.54632;
log_odds[8][2]=0.222332;
log_odds[8][3]=-0.01968;
log_odds[8][4]=0.233734;
log_odds[8][5]=-0.46956;
log_odds[8][6]=0.252154;
log_odds[8][7]=-0.18779;
log_odds[8][8]=0.266845;
log_odds[8][9]=-0.485;
log_odds[8][10]=-0.0275;
log_odds[8][11]=0.096405;
log_odds[8][12]=0.366055;*/
/* GLX */
log_odds[9][1]=0.0;
log_odds[9][2]=0.0;
log_odds[9][3]=0.0;
log_odds[9][4]=0.0;
log_odds[9][5]=0.0;
log_odds[9][6]=0.0;
log_odds[9][7]=0.0;
log_odds[9][8]=0.0;
log_odds[9][9]=0.0;
log_odds[9][10]=0.0;
log_odds[9][11]=0.0;
log_odds[9][12]=0.0;
/* GLY */
log_odds[10][1]=0.98;
log_odds[10][2]=-0.80;
log_odds[10][3]=0.90;
log_odds[10][4]=-1.18;
log_odds[10][5]=-0.49;
log_odds[10][6]=-1.52;
log_odds[10][7]=-0.98;
log_odds[10][8]=-1.76;
/*log_odds[10][1]=0.113474;
log_odds[10][2]=0.855742;
log_odds[10][3]=0.279528;
log_odds[10][4]=0.867145;
log_odds[10][5]=-0.41631;
log_odds[10][6]=0.142853;
log_odds[10][7]=-0.04929;
log_odds[10][8]=0.415072;
log_odds[10][9]=-1.37402;
log_odds[10][10]=-1.02188;
log_odds[10][11]=-0.95246;
log_odds[10][12]=-0.80437;*/
/* HIS */
log_odds[11][1]=0.10;
log_odds[11][2]=-0.00;
log_odds[11][3]=-0.18;
log_odds[11][4]=0.10;
log_odds[11][5]=0.25;
log_odds[11][6]=0.12;
log_odds[11][7]=-0.61;
log_odds[11][8]=-0.18;
/*log_odds[11][1]=-0.22366;
log_odds[11][2]=0.264696;
log_odds[11][3]=-0.04185;
log_odds[11][4]=0.345092;
log_odds[11][5]=-0.34513;
log_odds[11][6]=-0.62443;
log_odds[11][7]=0.134876;
log_odds[11][8]=0.021528;
log_odds[11][9]=0.311136;
log_odds[11][10]=-0.0905;
log_odds[11][11]=0.241285;
log_odds[11][12]=-0.44836;*/
/* ILE */
log_odds[12][1]=-0.65;
log_odds[12][2]=0.96;
log_odds[12][3]=-1.62;
log_odds[12][4]=-1.14;
log_odds[12][5]=-0.10;
log_odds[12][6]=0.67;
log_odds[12][7]=-0.82;
log_odds[12][8]=0.13;
/*log_odds[12][1]=0.299592;
log_odds[12][2]=-0.31067;
log_odds[12][3]=-1.31037;
log_odds[12][4]=-0.41259;
log_odds[12][5]=0.524639;
log_odds[12][6]=0.311658;
log_odds[12][7]=-0.62281;
log_odds[12][8]=-0.9593;
log_odds[12][9]=0.522446;
log_odds[12][10]=0.538113;
log_odds[12][11]=-0.5319;
log_odds[12][12]=-0.73604;*/
/* LEU */
log_odds[13][1]=-0.22;
log_odds[13][2]=0.20;
log_odds[13][3]=-0.72;
log_odds[13][4]=-0.46;
log_odds[13][5]=0.16;
log_odds[13][6]=0.31;
log_odds[13][7]=-0.08;
log_odds[13][8]=0.38;
/*log_odds[13][1]=-0.12288;
log_odds[13][2]=-0.03554;
log_odds[13][3]=-0.38291;
log_odds[13][4]=-1.09458;
log_odds[13][5]=0.339093;
log_odds[13][6]=0.333791;
log_odds[13][7]=-0.3885;
log_odds[13][8]=-0.56069;
log_odds[13][9]=0.549114;
log_odds[13][10]=0.174583;
log_odds[13][11]=-0.04976;
log_odds[13][12]=-0.31249;*/
/* LYS */
log_odds[14][1]=0.13;
log_odds[14][2]=-0.05;
log_odds[14][3]=0.18;
log_odds[14][4]=0.04;
log_odds[14][5]=0.07;
log_odds[14][6]=-0.15;
log_odds[14][7]=-0.25;
log_odds[14][8]=-0.18;
/*log_odds[14][1]=0.07382;
log_odds[14][2]=0.124959;
log_odds[14][3]=0.255622;
log_odds[14][4]=-0.07875;
log_odds[14][5]=-0.22998;
log_odds[14][6]=0.20367;
log_odds[14][7]=-0.01567;
log_odds[14][8]=0.218921;
log_odds[14][9]=-0.37877;
log_odds[14][10]=-0.00924;
log_odds[14][11]=-0.08028;
log_odds[14][12]=0.072263;*/
/* MET */
log_odds[15][1]=-0.19;
log_odds[15][2]=-0.04;
log_odds[15][3]=-0.18;
log_odds[15][4]=0.07;
log_odds[15][5]=0.24;
log_odds[15][6]=-0.07;
log_odds[15][7]=0.10;
log_odds[15][8]=-0.01;
/*log_odds[15][1]=0.243481;
log_odds[15][2]=-0.56745;
log_odds[15][3]=0.512294;
log_odds[15][4]=-0.3329;
log_odds[15][5]=-0.10944;
log_odds[15][6]=0.335178;
log_odds[15][7]=-0.18644;
log_odds[15][8]=-0.45394;
log_odds[15][9]=-0.02202;
log_odds[15][10]=0.543695;
log_odds[15][11]=-0.26544;
log_odds[15][12]=0.259942;*/
/* PHE */
log_odds[16][1]=-0.29;
log_odds[16][2]=-0.55;
log_odds[16][3]=-0.11;
log_odds[16][4]=0.96;
log_odds[16][5]=-0.09;
log_odds[16][6]=-0.59;
log_odds[16][7]=0.42;
log_odds[16][8]=-0.70;
/*log_odds[16][1]=-0.40378;
log_odds[16][2]=0.020032;
log_odds[16][3]=0.0;
log_odds[16][4]=-0.81586;
log_odds[16][5]=0.190362;
log_odds[16][6]=0.517195;
log_odds[16][7]=-1.71923;
log_odds[16][8]=-1.98673;
log_odds[16][9]=0.881316;
log_odds[16][10]=1.017232;
log_odds[16][11]=-0.48889;
log_odds[16][12]=-0.40534;*/
/* PRO */
log_odds[17][1]=-0.12;
log_odds[17][2]=-0.12;
log_odds[17][3]=-0.56;
log_odds[17][4]=-0.16;
log_odds[17][5]=0.14;
log_odds[17][6]=0.49;
log_odds[17][7]=-0.27;
log_odds[17][8]=0.27;
/*log_odds[17][1]=1.042648;
log_odds[17][2]=0.327028;
log_odds[17][3]=0.221148;
log_odds[17][4]=-0.0552;
log_odds[17][5]=0.839108;
log_odds[17][6]=-0.14926;
log_odds[17][7]=-0.00759;
log_odds[17][8]=-0.80118;
log_odds[17][9]=-0.87278;
log_odds[17][10]=-2.63085;
log_odds[17][11]=-1.12635;
log_odds[17][12]=-1.8901;*/
/* SER */
log_odds[18][1]=-1.36;
log_odds[18][2]=-1.28;
log_odds[18][3]=-0.49;
log_odds[18][4]=-0.51;
log_odds[18][5]=0.05;
log_odds[18][6]=0.03;
log_odds[18][7]=0.89;
log_odds[18][8]=0.86;
/*log_odds[18][1]=-0.14292;
log_odds[18][2]=-0.32174;
log_odds[18][3]=0.064857;
log_odds[18][4]=0.046338;
log_odds[18][5]=-0.1014;
log_odds[18][6]=-0.74087;
log_odds[18][7]=0.548317;
log_odds[18][8]=0.048196;
log_odds[18][9]=-0.16003;
log_odds[18][10]=-0.75697;
log_odds[18][11]=0.388581;
log_odds[18][12]=0.055235;*/
/* THR */
log_odds[19][1]=0.25;
log_odds[19][2]=0.53;
log_odds[19][3]=-0.23;
log_odds[19][4]=-0.11;
log_odds[19][5]=0.09;
log_odds[19][6]=-0.05;
log_odds[19][7]=-0.70;
log_odds[19][8]=-0.57;
/*log_odds[19][1]=0.113274;
log_odds[19][2]=0.269493;
log_odds[19][3]=0.129996;
log_odds[19][4]=0.319116;
log_odds[19][5]=0.034358;
log_odds[19][6]=0.342841;
log_odds[19][7]=0.254082;
log_odds[19][8]=0.089238;
log_odds[19][9]=-0.48436;
log_odds[19][10]=-0.18101;
log_odds[19][11]=-0.69141;
log_odds[19][12]=-0.60786;*/
/* TRP */
log_odds[20][1]=0.09;
log_odds[20][2]=-0.39;
log_odds[20][3]=-0.03;
log_odds[20][4]=0.78;
log_odds[20][5]=0.03;
log_odds[20][6]=-0.52;
log_odds[20][7]=0.15;
log_odds[20][8]=-1.06;
/*log_odds[20][1]=-1.5419;
log_odds[20][2]=0.268211;
log_odds[20][3]=0.0;
log_odds[20][4]=0.125463;
log_odds[20][5]=-0.3236;
log_odds[20][6]=-0.15092;
log_odds[20][7]=0.0;
log_odds[20][8]=-1.0454;
log_odds[20][9]=0.516391;
log_odds[20][10]=1.635786;
log_odds[20][11]=-0.05839;
log_odds[20][12]=-0.3803;*/
/* TYR */
log_odds[21][1]=-0.59;
log_odds[21][2]=-0.79;
log_odds[21][3]=-0.58;
log_odds[21][4]=1.48;
log_odds[21][5]=-0.52;
log_odds[21][6]=-0.51;
log_odds[21][7]=0.26;
log_odds[21][8]=-1.20;
/*log_odds[21][1]=-0.98667;
log_odds[21][2]=-0.38054;
log_odds[21][3]=-1.7857;
log_odds[21][4]=-2.49737;
log_odds[21][5]=-0.26737;
log_odds[21][6]=0.404307;
log_odds[21][7]=-1.60897;
log_odds[21][8]=-1.02917;
log_odds[21][9]=0.809251;
log_odds[21][10]=1.16139;
log_odds[21][11]=-0.28332;
log_odds[21][12]=0.537826;*/
/* VAL */
log_odds[22][1]=0.21;
log_odds[22][2]=0.80;
log_odds[22][3]=-0.87;
log_odds[22][4]=-1.16;
log_odds[22][5]=0.02;
log_odds[22][6]=0.31;
log_odds[22][7]=-1.02;
log_odds[22][8]=-0.15;
/*log_odds[22][1]=0.283564;
log_odds[22][2]=-0.55722;
log_odds[22][3]=-0.17062;
log_odds[22][4]=-0.37146;
log_odds[22][5]=0.460409;
log_odds[22][6]=0.035254;
log_odds[22][7]=-0.05105;
log_odds[22][8]=-0.62603;
log_odds[22][9]=0.513222;
log_odds[22][10]=-0.54469;
log_odds[22][11]=0.004303;
log_odds[22][12]=-0.69491;*/
/* IVA */
log_odds[23][1]=0.0;
log_odds[23][2]=0.0;
log_odds[23][3]=0.0;
log_odds[23][4]=0.0;
log_odds[23][5]=0.0;
log_odds[23][6]=0.0;
log_odds[23][7]=0.0;
log_odds[23][8]=0.0;
log_odds[23][9]=0.0;
log_odds[23][10]=0.0;
log_odds[23][11]=0.0;
log_odds[23][12]=0.0;
/* STA */
log_odds[24][1]=0.0;
log_odds[24][2]=0.0;
log_odds[24][3]=0.0;
log_odds[24][4]=0.0;
log_odds[24][5]=0.0;
log_odds[24][6]=0.0;
log_odds[24][7]=0.0;
log_odds[24][8]=0.0;
log_odds[24][9]=0.0;
log_odds[24][10]=0.0;
log_odds[24][11]=0.0;
log_odds[24][12]=0.0;
/* ACE */
log_odds[25][1]=0.0;
log_odds[25][2]=0.0;
log_odds[25][3]=0.0;
log_odds[25][4]=0.0;
log_odds[25][5]=0.0;
log_odds[25][6]=0.0;
log_odds[25][7]=0.0;
log_odds[25][8]=0.0;
log_odds[25][9]=0.0;
log_odds[25][10]=0.0;
log_odds[25][11]=0.0;
log_odds[25][12]=0.0;

/* THERMODYNAMIC PARAMETERS INITIALIZATION ------------------------------- */
aCp=0.44;
bCp=-0.26;
adeltaH=-8.44;
bdeltaH=31.4;
OHCp=0.17;
TsPolar=335.15;
TsApolar=385.15;
dSbb_length_corr=-0.12;
exposed_criteria=0.1;
W_Sconf=1.0;
Current_Temp = 298.15;
ASA_exposed_amide=7.5;
/* END THERMODYNAMIC PARAMETERS INITIALIZATION ------------------------------ */

}
