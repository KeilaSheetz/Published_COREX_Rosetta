/**** Ensemble Generator Application ****/
/**** Changed int2bin ****/

/*** STW 8/2002, RSS ***/

#define         Str_Len             100
#define         Default_FName       "PDB Filename"
#define         Default_Answer      "No"
#define         Def_Slice_Width     0.25
#define         Def_Solv_Rad        1.40
#define         Def_window_size     5
#define         Def_Min_Res_Win     4
#define         Max_Res             500
#define         Max_Atoms           3000
#define         Num_AA	            30
#define         Max_Nunits          100                         /* This also is the max number of windows for a partition */
#define         Max_WindowSize      30                         /* This also is the max number of partitions, Npartition */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <stdlib.h>
#include "asa_types.h"
typedef char bstr[Str_Len];

/* Global symbols */


float *x_coord;                                        /* main */
float *y_coord;
float *z_coord;
float *Atom_Area;
float *Atom_Area_Native;
bstr instring;
int window_size;
int Min_Res_Win;
int partition;
bstr out_fnameN;
FILE *outfile2;
char Atom_type[Max_Atoms][6];

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

FILE *infile;                                          /* Open_File */
bstr inFName;

int Total_Res;                                       /* read_pdb_file */
int OTnum;
char resname[Max_Res][8];
char atomname[Max_Atoms][6];
int Begin_Atom_Res[Max_Res];
int Last_Atom_Res[Max_Res];
int First_Res_Num;
int Total_Atoms;
int Last_Res_Num;
int pdb_resnum[Max_Res];
char pdb_molecule[Max_Res][6];

int Npartition;                                      /* Gen_Partition */
int First_Residue[Max_WindowSize][Max_Nunits];
int Last_Residue[Max_WindowSize][Max_Nunits];
int Nunits[Max_WindowSize];
long N_states[Max_WindowSize];

bstr out_fname;                                        /* Write_Atom_Info */
FILE *outfile1;

float ASA_N_Apolar;                                    /* Native_State */
float ASA_N_Polar;
float ASA_N_Apolar_unit[Max_Nunits];
float ASA_N_Polar_unit[Max_Nunits];
float ASA_U_Apolar_unit[Max_Nunits];
float ASA_U_Polar_unit[Max_Nunits];
float Fraction_exposed_Native[Max_Res];

char stateflag[Max_Nunits];                            /* int2bin */

int num_residues;                                      /* load_atoms_range */

float Sconf;                                           /* calc_ASA_dSconf */
float Fraction_amide_exposed[Max_Res];
float delASA_ap;
float delASA_pol;

/* Subroutines */

extern void Ecsact_Par (void);
extern void Open_File (void);
extern void read_pdb_file (asa_atom_rec_p *atom_list_p);
extern void Gen_Partition (void);
void Write_Atom_Info(int num_atoms, char label_nam[]);
extern void Native_State (void);
extern void int2bin (long k);
extern void calc_ASA_dSconf (int num_atoms);

void ZotNewLine();
int load_atoms(asa_atom_rec_p *atom_list_p);
int load_atoms_range(asa_atom_rec_p *atom_list_p);

/************************/

int main (int argc,char**argv)

{
asa_atom_rec_p atom_list;
int npointsmax;
int num_atoms;
char Answer[10];
int end;
double slice_width;
double solvent_radius;
int i;
char label_nam[7];
long Final_State;
long j;
float Fraction_Folded;
int jj;

npointsmax=Max_Atoms;
x_coord=(float*)malloc(npointsmax*sizeof(float));
y_coord=(float*)malloc(npointsmax*sizeof(float));
z_coord=(float*)malloc(npointsmax*sizeof(float));
Atom_Area=(float*)malloc(npointsmax*sizeof(float));
Atom_Area_Native=(float*)malloc(npointsmax*sizeof(float));

Ecsact_Par();
Open_File(); 
read_pdb_file(&atom_list);
fclose(infile);
num_atoms=load_atoms(&atom_list);
printf("Processing %d atoms\n",num_atoms);
if (!asa_assign_radii(atom_list,num_atoms,NULL,0))
   {
    fprintf(stdout,"Radii assignment failure.\n");
    exit(1);
   }

printf("\nEnter Number of Residues per Folding Unit [%d] > ",Def_window_size); 
    fgets(instring,Str_Len,stdin);
    ZotNewLine(instring);
    if (!strlen(instring))
        window_size=Def_window_size;
    else
        sscanf(instring,"%d",&window_size);

printf("\nEnter Minimum Number of Residues per Folding Unit [%d] > ",Def_Min_Res_Win); 
    fgets(instring,Str_Len,stdin);
    ZotNewLine(instring);
    if (!strlen(instring))
        Min_Res_Win=Def_Min_Res_Win;
    else
       sscanf(instring,"%d",&Min_Res_Win);

Gen_Partition();

end=1;
while(end)
     {
      printf("\nRun time too long? [%s] > ",Default_Answer); 
      fgets(instring,Str_Len,stdin);
      ZotNewLine(instring);
      if (!strlen(instring))
          strcpy(Answer,Default_Answer);
      else
          strcpy(Answer,instring);
      if ((strncmp(Answer,"Y",1)==0)||(strncmp(Answer,"y",1)==0))
         {
          printf("\nEnter Number of Residues per Folding Unit [%d] > ",Def_window_size); 
          fgets(instring,Str_Len,stdin);
          ZotNewLine(instring);
          if (!strlen(instring))
              window_size=Def_window_size;
          else
              sscanf(instring,"%d",&window_size);

          printf("\nEnter Minimum Number of Residues per Folding Unit [%d] > ",Def_Min_Res_Win); 
          fgets(instring,Str_Len,stdin);
          ZotNewLine(instring);
          if (!strlen(instring))
              Min_Res_Win=Def_Min_Res_Win;
          else
              sscanf(instring,"%d",&Min_Res_Win);

          Gen_Partition();
         }
      else
         {
          end=0;
         }
     }

solvent_radius=Def_Solv_Rad;
slice_width=Def_Slice_Width;
printf("\nrunning...");
fflush(stdout);
asa_calculate(atom_list,num_atoms,slice_width,solvent_radius);
for(i=0;i<num_atoms;i++) Atom_Area_Native[i]= atom_list[i].area;
Write_Atom_Info(num_atoms,label_nam);
sprintf(out_fname,"%s.%d.%d",inFName,window_size,Min_Res_Win);
printf("\n\nOut file = %s",out_fname);    		
outfile2=fopen(out_fname,"w");
for(partition=1;partition<=Npartition;partition++)
   {
    printf("\nNumber of intermediate states in partition %d = %ld\n",partition,N_states[partition]-2);
    printf("                                      ");
    Native_State();
    Final_State=N_states[partition]-1;
    if (partition==Npartition) Final_State=N_states[partition];
    for(j=1;j<Final_State;j++)
       {
        int2bin(j);
        printf("\r");
        printf("State %7ld	Code: %16s  ",j,stateflag);
        num_atoms=load_atoms_range(&atom_list);
        if (!asa_assign_radii(atom_list,num_atoms,NULL,0))
           {
            fprintf(stdout,"Radii assignment failure.\n");
            exit(1);
           }
        Fraction_Folded=(float)(num_residues)/(float)(Total_Res); 
        solvent_radius=Def_Solv_Rad;
        slice_width=Def_Slice_Width;
        fflush(stdout);
        asa_calculate(atom_list,num_atoms,slice_width,solvent_radius);
        for(i=0;i<num_atoms;i++)
           {
            Atom_Area[i]= atom_list[i].area;
            strncpy(Atom_type[i],atom_list[i].type,4);
           }
        calc_ASA_dSconf(num_atoms);
        fprintf(outfile2,"%d %ld %ld%8.4f%9.2f%9.2f%9.2f",partition,j,j,Fraction_Folded,Sconf,delASA_ap,delASA_pol);
        for(i=1;i<=Nunits[partition];i++)       
           {
            if (stateflag[i-1]=='0')
               {
                for(jj=First_Residue[partition][i];jj<=Last_Residue[partition][i];jj++) fprintf(outfile2,"%8.3f",Fraction_amide_exposed[jj]);
               }  
           }
        fprintf(outfile2,"\n");
       }
   }
fclose(outfile2);

asa_destroy_atom_list(atom_list);
printf("\n");
return(0);

}

/*********************Open_File ******************/

void Open_File(void)
{
printf("\nEnter pdb file name [%s] > ",Default_FName);
fgets(inFName,Str_Len,stdin);
ZotNewLine(inFName);
if (!strlen(inFName)) strcpy(inFName,Default_FName);
if ((infile=fopen(inFName, "r"))==NULL) 
   {
    fprintf(stderr,"\nfopen() error #%d on file %s\n",errno,inFName);
    exit(1);
   }
}

/****************** read_pdb_file ***************/

void read_pdb_file (asa_atom_rec_p *atom_list_p)
{
int end;
int atomnum;
int resnum,oldresnum;
int resnum_previous,resnum_old;
char molecule[2];
char string[120];
char junk[10],junk2[10];

end=1;
Total_Res=0;
atomnum=0;
oldresnum=-100;
resnum=0;
resnum_previous=10000;
OTnum=0;
while(end)
     {
      fgets(string,100,infile);
      if(feof(infile))
	break;
      if (strncmp(string,"END",3) == 0)
         {
          break;
         }
      if (strncmp(string,"ATOM",4) == 0)
         {
          if (string[13] == 'H') continue;                                                             /* stop-gap hydrogen eliminator */
          sscanf(&string[12],"%3s", junk2);
          sscanf(&string[22],"%4s", junk);
          sscanf(&string[21],"%1s", molecule);
          resnum_old=atoi(junk);
          if (resnum_previous!=resnum_old)
             {
              resnum_previous=resnum_old;
              resnum++;
             }
          pdb_resnum[resnum]=resnum_old;
          sscanf(&string[17],"%3s", resname[resnum]);
          sscanf(&string[30],"%f", &x_coord[atomnum]);
          sscanf(&string[38],"%f", &y_coord[atomnum]);
          sscanf(&string[46],"%f", &z_coord[atomnum]);
          strcpy(pdb_molecule[resnum]," ");
          if (string[21] == ' ')
             {
              strcat(pdb_molecule[resnum]," ");
             }
          else
             {
              sscanf(&string[21],"%1s", molecule);
              strcat(pdb_molecule[resnum],molecule);
             }
          strcat(pdb_molecule[resnum]," ");
          strcpy(atomname[atomnum]," ");
          strcat(atomname[atomnum],junk2);
          strcat(atomname[atomnum],"  ");
          if (resnum!=oldresnum)
             {
              Begin_Atom_Res[resnum]=atomnum;
              Last_Atom_Res[resnum-1]=atomnum-1;
             }
          if (strncmp(atomname[atomnum]," OXT",4)==0||strncmp(atomname[atomnum]," OT",3)==0) OTnum++;
          if (strncmp(atomname[atomnum]," CA",3)==0) Total_Res++;
          atomnum++;
          oldresnum=resnum;
          if (atomnum==1) First_Res_Num=resnum;
         }
     }
if (Total_Res<resnum)
   {
    printf("\nAbnormality in input. Missing alpha carbon (CA) in pdb file.\nCheck pdb file. PROGRAM WILL CONTINUE.\n");
    printf("\nPress <return> to continue "); 
    fgets(instring, Str_Len, stdin);
    ZotNewLine(instring);
    printf("\n");
    Total_Res=resnum;
   }
if (Total_Res>resnum)
   {
    printf("\nAbnormality in input. Two sequentially listed residues in pdb file have\nsame residue number. This likely will affect residue specific calculations.\nCheck pdb file. PROGRAM WILL CONTINUE.\n");
    printf("\nPress <return> to continue "); 
    fgets(instring, Str_Len, stdin);
    ZotNewLine(instring);
    printf("\n");
   }
Total_Atoms=atomnum;
Last_Atom_Res[resnum]=atomnum-1;
Last_Res_Num=First_Res_Num+Total_Res-1;
printf("Total Number of residues = %d. Last atom %s %d\n",Total_Res,atomname[atomnum-1],atomnum-1);
printf("First Residue %s %d %s, Last Residue %s %d %s\n",resname[First_Res_Num],pdb_resnum[First_Res_Num],pdb_molecule[First_Res_Num],resname[Last_Res_Num],pdb_resnum[Last_Res_Num],pdb_molecule[Last_Res_Num]);

}

/************************** load_atoms **********************/

int load_atoms(asa_atom_rec_p *atom_list_p)
{
int num_atoms,i;
asa_atom_rec_p atom_list;

num_atoms=Total_Atoms;
asa_create_atom_list(atom_list_p,num_atoms); 
atom_list = *atom_list_p;
for(i=0;i<num_atoms;i++)
   {
    strncpy(atom_list[i].type,atomname[i],4); 
    atom_list[i].x = (double) x_coord[i];
    atom_list[i].y = (double) y_coord[i];
    atom_list[i].z = (double) z_coord[i]; 
    atom_list[i].area=0.0; 
   }
return num_atoms;

}

/************************** Gen_Partition ******************/

void Gen_Partition (void)
{
int partition_index;
int end;
int i;
int k;
long sum;
float time_hours,time_days,time_minutes;

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
    if (Last_Residue[partition_index][1] - First_Residue[partition_index][1]< Min_Res_Win-1)
       {
        Last_Residue[partition_index][1] = Last_Residue[partition_index][2];
        for(i=2;i< Nunits[partition_index];i++)
           {
            First_Residue[partition_index][i] = First_Residue[partition_index][i+1];
            Last_Residue[partition_index][i] = Last_Residue[partition_index][i+1];
           }
        Nunits[partition_index] = Nunits[partition_index] - 1;
       }
    if (Last_Residue[partition_index][Nunits[partition_index]] - First_Residue[partition_index][Nunits[partition_index]]< Min_Res_Win-1)
       {
        Last_Residue[partition_index][Nunits[partition_index]-1] = Last_Residue[partition_index][Nunits[partition_index]];
        Nunits[partition_index] = Nunits[partition_index] - 1;
       }
   }

sum=0;
printf("\nTotal Number of Partitions = %d\n",Npartition);
for(k=1;k<=Npartition;k++)
   {
    N_states[k]=(pow((double)2.0,(double)Nunits[k]));
    printf("There are %ld intermediate states in partition %d\n", N_states[k]-2,k);
    sum=sum+N_states[k]-2;
   }

printf("Total Number of States in Ensemble = %ld\n",sum+2);
time_hours=(sum+2)*72.0/1000000.0;
time_days=time_hours/24.0;
time_minutes=time_hours*60.0;
printf("Approximate run time = %f minutes = %f hours = %f days\n",time_minutes,time_hours,time_days);

}

/***************Write_Atom_Info*****************/

void Write_Atom_Info(int num_atoms, char label_nam[])
{
int i,j;

strcpy(out_fname,inFName);
strcat(out_fname,".info");
outfile1=fopen(out_fname,"w");
fprintf(outfile1,"%d\n",num_atoms);
fprintf(outfile1,"ResName ResNum AtomNum AtomName Nat.Area pdbResNum\n");
for(j=First_Res_Num;j<= Last_Res_Num;j++)
   {
    for(i=Begin_Atom_Res[j];i<=Last_Atom_Res[j];i++)
       {
        strncpy(label_nam,atomname[i],4);
        label_nam[4]='\0';
        fprintf(outfile1,"%3s%5d%5d%4s%6.2f%6d%s\n",resname[j],j,i,label_nam,Atom_Area_Native[i],pdb_resnum[j],pdb_molecule[j]);
       }
   }

fclose(outfile1);

}

/****************************** Native_State ****************************/

void Native_State ()
{
int i,k;
int j,jj;
char junk[7];
float ASA_side_chain;

ASA_N_Apolar=0.0;
ASA_N_Polar=0.0;
for(i=1;i<=Nunits[partition];i++)
   {
    ASA_N_Apolar_unit[i]=0.0;
    ASA_N_Polar_unit[i]=0.0;
    ASA_U_Apolar_unit[i]=0.0;
    ASA_U_Polar_unit[i]=0.0;
   }

strcpy(junk,"0000");
for(i=1;i<=Nunits[partition];i++)
   {
    for(j=First_Residue[partition][i];j<=Last_Residue[partition][i];j++)
       {
        ASA_side_chain=0.0;
        for(k=Begin_Atom_Res[j];k<=Last_Atom_Res[j];k++)
           {
            strncpy(junk,atomname[k],4);
            if (junk[1]=='C')
               {				
                ASA_N_Apolar_unit[i]+=Atom_Area_Native[k];
                ASA_N_Apolar=ASA_N_Apolar+Atom_Area_Native[k];
               }
            else
               {
                ASA_N_Polar_unit[i]+=Atom_Area_Native[k];
                ASA_N_Polar=ASA_N_Polar+Atom_Area_Native[k];
               }
            if (strncmp(junk," N ",3)&&strncmp(junk," CA",3)&&strncmp(junk," C ",3)&&strncmp(junk," O ",3)) ASA_side_chain=ASA_side_chain+Atom_Area_Native[k];
           }
        for(jj=1;jj<=NAATotal;jj++)
           {
            if (strncmp(resname[j],aa3[jj],3)==0)
               {
                ASA_U_Apolar_unit[i]=ASA_U_Apolar_unit[i]+ASAexapol[jj];
                ASA_U_Polar_unit[i]=ASA_U_Polar_unit[i]+ASAexpol[jj];
                Fraction_exposed_Native[j]=ASA_side_chain/ASAsc[jj];
               }
           }
       }
   }
ASA_U_Polar_unit[1]=ASA_U_Polar_unit[1]+45.0;
j=Nunits[partition];
ASA_U_Apolar_unit[j]=ASA_U_Apolar_unit[j]+30.0;
ASA_U_Polar_unit[j]=ASA_U_Polar_unit[j]+30.0*OTnum;

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

/********************************** load_atoms_range ****************************/

int load_atoms_range(asa_atom_rec_p *atom_list_p)
{
int num_atoms,i,j,k,atom_index[Max_Atoms];
asa_atom_rec_p atom_list;
char junk[7],junk1[7];
strcpy(junk,"0000");
strcpy(junk1,"0000");

num_atoms=0;
num_residues=0;
for(i=1;i<=Nunits[partition];i++)
   {
    if (stateflag[i-1] == '0')
       {
        for(j=First_Residue[partition][i];j<=Last_Residue[partition][i];j++)
           {
            for(k=Begin_Atom_Res[j];k<=Last_Atom_Res[j];k++)
               {	
                num_atoms++;
               }
           }
       }
   }

num_atoms=num_atoms+1;
atom_list = *atom_list_p; 
num_atoms=0;
for(i=1;i<=Nunits[partition];i++)
   {
    if (stateflag[i-1] == '0')
       {
        for(j=First_Residue[partition][i];j<=Last_Residue[partition][i];j++)
           {
            num_residues++;
            for(k=Begin_Atom_Res[j];k<=Last_Atom_Res[j];k++)
               {				
                strncpy(atom_list[num_atoms].type, atomname[k], 4);  
                atom_list[num_atoms].x = (double) x_coord[k];
                atom_list[num_atoms].y = (double) y_coord[k];
                atom_list[num_atoms].z = (double) z_coord[k]; 
                atom_list[num_atoms].area=0.0;
                atom_index[num_atoms]=k;
                num_atoms++;
               }
           }
       }
   }
num_atoms=num_atoms; 
strncpy(junk,atom_list[0].type,4);
k=atom_index[num_atoms-1];
strncpy(junk1,atom_list[num_atoms-1].type,4);
return num_atoms;

}

/******************************** calc_ASA_dSconf *********************************/

void calc_ASA_dSconf (int num_atoms)
{
int i,j,k,atom_number,jj;
char junk[7];
float Sum_U_Apolar,Sum_U_Polar,ASA_side_chain,ASA_Amide;
float ASA_State_Apolar,ASA_State_Polar;
float Fraction_exposed[Max_Res];

ASA_State_Apolar=0.0;
ASA_State_Polar=0.0;
Sum_U_Apolar=0.0;
Sum_U_Polar=0.0;
 
for(i=0;i<num_atoms;i++)
   { 			                                 
    strncpy(junk,Atom_type[i],4);
    if (junk[1]=='C')
       {
        ASA_State_Apolar=ASA_State_Apolar+Atom_Area[i];
       }
    else
       {
        ASA_State_Polar=ASA_State_Polar+Atom_Area[i];
       } 
   }

Sconf=0.0;
for(i=1;i<=Nunits[partition];i++)
   {
    if (stateflag[i-1]=='1')
       {
        Sum_U_Apolar=Sum_U_Apolar+ASA_U_Apolar_unit[i];
        Sum_U_Polar=Sum_U_Polar+ASA_U_Polar_unit[i];
        for(j=First_Residue[partition][i];j<=Last_Residue[partition][i];j++)
           {
            for(jj=1;jj<=NAATotal;jj++)
               {
                if (strncmp(resname[j],aa3[jj],3)==0)
                   {
                    Sconf=Sconf+dSexu[jj]+dSbb[jj]+dSbb_length_corr+(1.0-Fraction_exposed_Native[j])*dSbuex[jj];
                   }			
               }
           }
       }
   }

strcpy(junk,"0000");
atom_number=0;
for(i=1;i<=Nunits[partition];i++)
   {
    if (stateflag[i-1] == '0')
       {
        for(j=First_Residue[partition][i];j<=Last_Residue[partition][i];j++)
           {
            ASA_side_chain=0.0;
            ASA_Amide=0.0;
            for(k=Begin_Atom_Res[j];k<=Last_Atom_Res[j];k++)
               {
                strncpy(junk,Atom_type[atom_number], 4); 
                if (strncmp(junk," N ",3)==0) ASA_Amide=ASA_Amide+Atom_Area[atom_number];
                if (strncmp(junk," N ",3)&&strncmp(junk," CA",3)&&strncmp(junk," C ",3)&& strncmp(junk," O ",3)) ASA_side_chain=ASA_side_chain+Atom_Area[atom_number];
                atom_number++;
               }
            for(jj=1;jj<=NAATotal;jj++)
               {
                if (strncmp(resname[j],aa3[jj],3)==0)
                   {
                    Fraction_exposed[j]=ASA_side_chain/ASAsc[jj];
                    Fraction_amide_exposed[j]=ASA_Amide/ASA_exposed_amide;
                    Sconf=Sconf+(Fraction_exposed[j]-Fraction_exposed_Native[j])*dSbuex[jj];
                   }
               }
           }
       }
   }

delASA_ap=ASA_State_Apolar+Sum_U_Apolar-ASA_N_Apolar;
delASA_pol=ASA_State_Polar+Sum_U_Polar-ASA_N_Polar;
	
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
