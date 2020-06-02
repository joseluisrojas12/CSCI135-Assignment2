#include<iostream> 
#include<fstream>
#include<cstring>
#include<cctype>
#include<cmath>

#include<string>

using namespace std;


class DNASequence{

private:
string firstdna;
string aFile;
int N;
int jamino;
int c;
int j;
int jj;
int aminocount;
int p;//for aminofromFile
int m;//for aminofromFile
char n1,n2,n3;
char na1,na2,na3;
char aminoAcid;
char sentenceNucleotideNew[40];
char sentenceCodon[20];
char amino[];
char firstdnastring[40];
int counter;

public:

int getLength();
void clear();

void appendCodon(char n1, char n2, char n3);
void appendAminoAcid(char aminoAcid);
void appendAminoAcidsFromFile(string aFile);//*************************

string getNucleotideSequence();
string getAminoAcidSequence();

};


int DNASequence::getLength()
{
 jj=j;
 return jj;
}

void DNASequence::clear()
{
 n1='\0';
 n2='\0';
 n3='\0';
 
 c=0;
 j=0;

 strcpy(sentenceNucleotideNew,"");
}


void DNASequence::appendCodon(char n1, char n2, char n3)
{c=1;
 this->n1=n1;
 this->n2=n2;
 this->n3=n3;

 //if(c==0)
   //{
    //strcpy(sentenceNucleotideNew,"");
   //}
 
  if(((n1=='a'||n1=='A')||
      (n1=='c'||n1=='C')||
      (n1=='t'||n1=='T')||
      (n1=='g'||n1=='G'))&&
     ((n2=='a'||n2=='A')||
      (n2=='c'||n2=='C')||
      (n2=='t'||n2=='T')||
      (n2=='g'||n2=='G'))&&
     ((n3=='a'||n3=='A')||
      (n3=='c'||n3=='C')||
      (n3=='t'||n3=='T')||
      (n3=='g'||n3=='G')))
   
     
      {
       if(n1=='a')
        n1='A';
       if(n1=='c')
        n1='C';
       if(n1=='g')
        n1='G';
       if(n1=='t')
        n1='T';
       if(n2=='a')
        n2='A';
       if(n2=='c')
        n2='C';
       if(n2=='g')
        n2='G';
       if(n2=='t')
        n2='T';
       if(n3=='a')
        n3='A';
       if(n3=='c')
        n3='C';
       if(n3=='g')
        n3='G';
       if(n3=='t')
        n3='T';
       
       char word1[2]={""};
       char word2[2]={""};
       char word3[2]={""};

       word1[0]=n1;
       word2[0]=n2;
       word3[0]=n3;

       char sentenceNucleotide[40]={""};
        
       strcat(sentenceNucleotide,word1);
       strcat(sentenceNucleotide,word2);
       strcat(sentenceNucleotide,word3);
               
       strcat(sentenceNucleotideNew,sentenceNucleotide);

       j=strlen(sentenceNucleotideNew);//This is to get the lenght.

       c=1;

      }
  else
      {
       cout<<"";//The chain of nucleotides is not legal"<<endl;  Used for testing. 
       
      }


}


void DNASequence::appendAminoAcid(char aminoAcid)
{
  if(aminoAcid=='f'||aminoAcid=='F')
   {
    n1='T';
    n2='T';//char amino[]={"TTT"};
    n3='T';
   }
  else if(aminoAcid=='l'||aminoAcid=='L')
   {
    n1='T';
    n2='T';
    n3='A';//char amino[]={"TTA"};
   }
  else if(aminoAcid=='i'||aminoAcid=='I')
   {
    n1='A';
    n2='T';//char amino[]={"ATT"};
    n3='T';
   }
  else if(aminoAcid=='m'||aminoAcid=='M')//OJO:Este es START
   {
    n1='A';
    n2='T';//char amino[]={"ATG"};
    n3='G';
   }
  else if(aminoAcid=='v'||aminoAcid=='V')
   {
    n1='G';
    n2='T';//char amino[]={"GTT"};
    n3='T';
   }
  else if(aminoAcid=='s'||aminoAcid=='S')
   {
    n1='T';
    n2='C';//char amino[]={"TCT"};
    n3='T';
   }
  else if(aminoAcid=='p'||aminoAcid=='P')
   {
    n1='C';
    n2='C';//char amino[]={"CCT"};
    n3='T';
   }
  else if(aminoAcid=='t'||aminoAcid=='T')
   {
    n1='A';
    n2='C';
    n3='T';//char amino[]={"ACT"};
   }
  else if(aminoAcid=='a'||aminoAcid=='A')
   {
    n1='G';
    n2='C';//char amino[]={"GCT"};
    n3='T';
   }
  else if(aminoAcid=='y'||aminoAcid=='Y')
   {
    n1='T';
    n2='A';//char amino[]={"TAT"};
    n3='T';
   }
  else if(aminoAcid=='!'||aminoAcid=='!')//STOP
   {
    char amino[]={"TAA"};
   }
  else if(aminoAcid=='h'||aminoAcid=='H')
   {
    n1='C';
    n2='A';//char amino[]={"CAT"};
    n3='T';
   }
  else if(aminoAcid=='q'||aminoAcid=='Q')
   {
    n1='C';
    n2='A';//char amino[]={"CAA"};
    n3='A';
   }
  else if(aminoAcid=='n'||aminoAcid=='N')
   {
    n1='A';
    n2='A';//char amino[]={"AAT"};
    n3='T';
   }
  else if(aminoAcid=='k'||aminoAcid=='K')
   {
    n1='A';
    n2='A';//char amino[]={"AAA"};
    n3='A';
   }
  else if(aminoAcid=='d'||aminoAcid=='D')
   {
    n1='G';
    n2='A';//char amino[]={"GAT"};
    n3='T';
   }
  else if(aminoAcid=='e'||aminoAcid=='E')
   {
    n1='G';
    n2='A';
    n3='A';//char amino[]={"GAA"};
   }
  else if(aminoAcid=='c'||aminoAcid=='C')
   {
    n1='T';
    n2='G';//char amino[]={"TGC"};
    n3='C';
   }
  else if(aminoAcid=='w'||aminoAcid=='W')
   {
    n1='T';
    n2='G';
    n3='G';//char amino[]={"TGG"};
   }
  else if(aminoAcid=='r'||aminoAcid=='R')
   {
    n1='C';
    n2='G';
    n3='T';//char amino[]={"CGT"};
   }
  else if(aminoAcid=='s'||aminoAcid=='S')
   {
    n1='A';
    n2='G';
    n3='T';//char amino[]={"AGT"};
   }
  else if(aminoAcid=='g'||aminoAcid=='G')
   {
    n1='G';
    n2='G';//char amino[]={"GGT"};
    n3='T';
   }
  else
   {
    cout<<"";//"Is not a valid amino char"<<endl;   Used for testing.
    n1='\0';  
    n2='\0';
    n3='\0';
   }
 //cout<<n1<<n2<<n3<<endl;//Display the valid amino   OJO: Puesto ahorita??????

   if(((n1=='a'||n1=='A')||
      (n1=='c'||n1=='C')||
      (n1=='t'||n1=='T')||
      (n1=='g'||n1=='G'))&&
     ((n2=='a'||n2=='A')||
      (n2=='c'||n2=='C')||
      (n2=='t'||n2=='T')||
      (n2=='g'||n2=='G'))&&
     ((n3=='a'||n3=='A')||
      (n3=='c'||n3=='C')||
      (n3=='t'||n3=='T')||
      (n3=='g'||n3=='G')))
   
     
      {
       if(n1=='a')
        n1='A';
       if(n1=='c')
        n1='C';
       if(n1=='g')
        n1='G';
       if(n1=='t')
        n1='T';
       if(n2=='a')
        n2='A';
       if(n2=='c')
        n2='C';
       if(n2=='g')
        n2='G';
       if(n2=='t')
        n2='T';
       if(n3=='a')
        n3='A';
       if(n3=='c')
        n3='C';
       if(n3=='g')
        n3='G';
       if(n3=='t')
        n3='T';
       
       char word1[2]={""};
       char word2[2]={""};
       char word3[2]={""};

       word1[0]=n1;
       word2[0]=n2;
       word3[0]=n3;

       char sentenceNucleotide[40]={""};
        
       strcat(sentenceNucleotide,word1);
       strcat(sentenceNucleotide,word2);
       strcat(sentenceNucleotide,word3);
               
       strcat(sentenceNucleotideNew,sentenceNucleotide);

       j=strlen(sentenceNucleotideNew);//This is to get the lenght.

       c=1;

      }
  else
      {
       cout<<"";//The chain of nucleotides is not legal"<<endl;  Used for testing. 
       
      }


}


void DNASequence::appendAminoAcidsFromFile(string aFile)//()**********
{
 ifstream inputStream;//fstream inputStream;***************
 inputStream.open("amino_human_iris.txt");//*****************
 if(inputStream.fail())
  {
   cerr<<"No existe el file amino_human_iris.txt"<<endl;
  }
 else
  {//File found and ready to process.
   const int LENGTH=20;
   jamino=0;
   aminocount=0;
   inputStream>>firstdna;
   //cout<<"First amino chain in amino1.txt is: "<<firstdna<<endl; Used to check what is in th file.
   for(int i=0;i<LENGTH;i++,jamino++)
    {
     if((firstdna[i]=='f'||firstdna[i]=='F')||
        (firstdna[i]=='l'||firstdna[i]=='L')||
        (firstdna[i]=='i'||firstdna[i]=='I')||   
        (firstdna[i]=='m'||firstdna[i]=='M')||//OJO:Este es START
        (firstdna[i]=='v'||firstdna[i]=='V')||
        (firstdna[i]=='s'||firstdna[i]=='S')||
        (firstdna[i]=='p'||firstdna[i]=='P')||
        (firstdna[i]=='t'||firstdna[i]=='T')||
        (firstdna[i]=='a'||firstdna[i]=='A')||
        (firstdna[i]=='y'||firstdna[i]=='Y')||
        (firstdna[i]=='!'||firstdna[i]=='!')||//OJO: STOP
        (firstdna[i]=='h'||firstdna[i]=='H')||
        (firstdna[i]=='q'||firstdna[i]=='Q')||
        (firstdna[i]=='n'||firstdna[i]=='N')||
        (firstdna[i]=='k'||firstdna[i]=='K')||
        (firstdna[i]=='d'||firstdna[i]=='D')||
        (firstdna[i]=='e'||firstdna[i]=='E')||
        (firstdna[i]=='c'||firstdna[i]=='C')||
        (firstdna[i]=='w'||firstdna[i]=='W')||
        (firstdna[i]=='r'||firstdna[i]=='R')||
        (firstdna[i]=='s'||firstdna[i]=='S')||
        (firstdna[i]=='g'||firstdna[i]=='G'))
        {
         cout<<"";//"El Amino es legal"<<endl;  Used for testing.
         cout<<"";//Character "<<firstdna[i]<<" valido y filtrado"<<endl;   Used for testing.
         firstdnastring[jamino]=firstdna[i];
         aminocount=aminocount+1;
        }
        else
        {
         jamino=i+1;
         if(firstdna[i]=='\0')
         {
                    
          break;
         }
        }
    }
   for(int pe=0;pe<=jamino;pe++)
    {
     if(firstdnastring[pe]!='\0')
       {
        //cout<<firstdnastring[pe];  Used for testing filtered file.
        counter=pe;
        if(firstdnastring[counter]=='f'||firstdnastring[counter]=='F')
        {
         na1='T';
         na2='T';//char amino[]={"TTT"};
         na3='T';
        }
        else if(firstdnastring[counter]=='l'||firstdnastring[counter]=='L')
        {
         na1='T';
         na2='T';
         na3='A';//char amino[]={"TTA"};
        }
        else if(firstdnastring[counter]=='i'||firstdnastring[counter]=='I')
        {
         na1='A';
         na2='T';//char amino[]={"ATT"};
         na3='T';
        }
        else if(firstdnastring[counter]=='m'||firstdnastring[counter]=='M')//OJO:Este es START
        {
         na1='A';
         na2='T';//char amino[]={"ATG"};
         na3='G';
        }
        else if(firstdnastring[counter]=='v'||firstdnastring[counter]=='V')
        {
         na1='G';
         na2='T';//char amino[]={"GTT"};
         na3='T';
        }
        else if(firstdnastring[counter]=='s'||firstdnastring[counter]=='S')
        {
         na1='T';
         na2='C';//char amino[]={"TCT"};
         na3='T';
        }
        else if(firstdnastring[counter]=='p'||firstdnastring[counter]=='P')
        {
         na1='C';
         na2='C';//char amino[]={"CCT"};
         na3='T';
        }
        else if(firstdnastring[counter]=='t'||firstdnastring[counter]=='T')
        {
         na1='A';
         na2='C';
         na3='T';//char amino[]={"ACT"};
        }
        else if(firstdnastring[counter]=='a'||firstdnastring[counter]=='A')
        {
         na1='G';
         na2='C';//char amino[]={"GCT"};
         na3='T';
        }
        else if(firstdnastring[counter]=='y'||firstdnastring[counter]=='Y')
        {
         na1='T';
         na2='A';//char amino[]={"TAT"};
         na3='T';
        }
        else if(firstdnastring[counter]=='h'||firstdnastring[counter]=='H')
        {
         na1='C';
         na2='A';//char amino[]={"CAT"};
         na3='T';
        }
        else if(firstdnastring[counter]=='q'||firstdnastring[counter]=='Q')
        {
         na1='C';
         na2='A';//char amino[]={"CAA"};
         na3='A';
        }
        else if(firstdnastring[counter]=='!'||firstdnastring[counter]=='!')
        {
         na1='T';
         na2='A';//char amino[]={"AAT"};
         na3='A';
        }
        else if(firstdnastring[counter]=='n'||firstdnastring[counter]=='N')
        {
         na1='A';
         na2='A';//char amino[]={"AAT"};
         na3='T';
        }
        else if(firstdnastring[counter]=='k'||firstdnastring[counter]=='K')
        {
         na1='A';
         na2='A';//char amino[]={"AAA"};
         na3='A';
        }
        else if(firstdnastring[counter]=='d'||firstdnastring[counter]=='D')
        {
         na1='G';
         na2='A';//char amino[]={"GAT"};
         na3='T';
        }
        else if(firstdnastring[counter]=='e'||firstdnastring[counter]=='E')
        {
         na1='G';
         na2='A';
         na3='A';//char amino[]={"GAA"};
        }
        else if(firstdnastring[counter]=='c'||firstdnastring[counter]=='C')
        {
         na1='T';
         na2='G';//char amino[]={"TGC"};
         na3='C';
        }
        else if(firstdnastring[counter]=='w'||firstdnastring[counter]=='W')
        {
         na1='T';
         na2='G';
         na3='G';//char amino[]={"TGG"};
        }
        else if(firstdnastring[counter]=='r'||firstdnastring[counter]=='R')
        {
         na1='C';
         na2='G';
         na3='T';//char amino[]={"CGT"};
        }
        else if(firstdnastring[counter]=='s'||firstdnastring[counter]=='S')
        {
         na1='A';
         na2='G';
         na3='T';//char amino[]={"AGT"};
        }
        else if(firstdnastring[counter]=='g'||firstdnastring[counter]=='G')
        {
         na1='G';
         na2='G';//char amino[]={"GGT"};
         na3='T';
        }
        else
        {
         cout<<"";//"Is not a valid amino char"<<endl;   Used for testing.
         break;
    //na1='\0';
    //na2='\0';
    //na3='\0';
        }
   
        char wordaminofile1[2]={""};
        char wordaminofile2[2]={""};
        char wordaminofile3[2]={""};

        wordaminofile1[0]=na1;
        wordaminofile2[0]=na2;
        wordaminofile3[0]=na3;

        char sentenceNucleotideFile[40]={""};

        strcat(sentenceNucleotideFile,wordaminofile1);
        strcat(sentenceNucleotideFile,wordaminofile2);
        strcat(sentenceNucleotideFile,wordaminofile3);
        strcat(sentenceNucleotideNew,sentenceNucleotideFile);
       }

    }


   //cout<<"longitud de firstdnastring es: "<<aminocount<<endl;   Used for testing.
   j=j+aminocount*3;
  }


 n1='\0';
 n2='\0';
 n3='\0';   
}


string DNASequence::getNucleotideSequence()
{

 
 return sentenceNucleotideNew;
}


string DNASequence::getAminoAcidSequence()
{
 char sentenceAminoFinalNew[40]={""};
 char sentenceAminoOutput[40]={""};
 strcpy(sentenceAminoFinalNew,sentenceNucleotideNew);//To avoid segment violation.
 for(int counternucleotidefinal=0;counternucleotidefinal<j;counternucleotidefinal++) 
  {
    if((sentenceAminoFinalNew[counternucleotidefinal]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='G'))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='*';
        counternucleotidefinal=counternucleotidefinal+2;
       }
    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='T'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='C')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='F';
        counternucleotidefinal=counternucleotidefinal+2;
       }
    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='A'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='G')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='L';
        counternucleotidefinal=counternucleotidefinal+2;
       }
    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='T'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='C')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='L';
        counternucleotidefinal=counternucleotidefinal+2;
       }
    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='A'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='G')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='L';
        counternucleotidefinal=counternucleotidefinal+2;
       }
    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='T'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='C'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='A')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='I';
        counternucleotidefinal=counternucleotidefinal+2;
       }
    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='T'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='C')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='V';
        counternucleotidefinal=counternucleotidefinal+2;
       }
    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='A'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='G')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='V';
        counternucleotidefinal=counternucleotidefinal+2;
       }
    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='T'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='C')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='S';
        counternucleotidefinal=counternucleotidefinal+2;
       }
    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='A'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='G')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='S';
        counternucleotidefinal=counternucleotidefinal+2;
       }

    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='T'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='C')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='P';
        counternucleotidefinal=counternucleotidefinal+2;
       }
    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='A'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='G')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='P';
        counternucleotidefinal=counternucleotidefinal+2;
       }
    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='T'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='C')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='T';
        counternucleotidefinal=counternucleotidefinal+2;
       }
    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='A'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='G')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='T';
        counternucleotidefinal=counternucleotidefinal+2;
       }

    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='T'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='C')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='A';
        counternucleotidefinal=counternucleotidefinal+2;
       }

    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='A'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='G')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='A';
        counternucleotidefinal=counternucleotidefinal+2;
       }

    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='T'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='C')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='Y';
        counternucleotidefinal=counternucleotidefinal+2;
       }

    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='A'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='G'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='A')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='!';
        counternucleotidefinal=counternucleotidefinal+2;
       }

    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='T'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='C')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='H';
        counternucleotidefinal=counternucleotidefinal+2;
       }

    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='A'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='G')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='Q';
        counternucleotidefinal=counternucleotidefinal+2;
       }

    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='T'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='C')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='N';
        counternucleotidefinal=counternucleotidefinal+2;
       }

    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='A'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='G')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='K';
        counternucleotidefinal=counternucleotidefinal+2;
       }

    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='T'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='C')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='D';
        counternucleotidefinal=counternucleotidefinal+2;
       }

    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='A'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='G')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='E';
        counternucleotidefinal=counternucleotidefinal+2;
       }

    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='C'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='T')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='C';
        counternucleotidefinal=counternucleotidefinal+2;
       }

    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='G'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='T')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='G')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='W';
        counternucleotidefinal=counternucleotidefinal+2;
       }

    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='T'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='C')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='R';
        counternucleotidefinal=counternucleotidefinal+2;
       }

    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='A'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='A')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='R';
        counternucleotidefinal=counternucleotidefinal+2;
       }

    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='G'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='C')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='G')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='R';
        counternucleotidefinal=counternucleotidefinal+2;
       }

    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='T'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='A')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='C')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='S';
        counternucleotidefinal=counternucleotidefinal+2;
       }

    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='T'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='C')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='G';
        counternucleotidefinal=counternucleotidefinal+2;
       }

    else if(((sentenceAminoFinalNew[counternucleotidefinal]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='A'))||
       ((sentenceAminoFinalNew[counternucleotidefinal]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+1]=='G')&&
       (sentenceAminoFinalNew[counternucleotidefinal+2]=='G')))
       {
        sentenceAminoOutput[counternucleotidefinal/3]='G';
        counternucleotidefinal=counternucleotidefinal+2;
       }


    else
       {
        if(sentenceAminoFinalNew[counternucleotidefinal]=='\0')
          {
           sentenceAminoOutput[counternucleotidefinal]=sentenceAminoFinalNew[counternucleotidefinal];
          }
         else
          {
           sentenceAminoOutput[counternucleotidefinal/3]='A';
           counternucleotidefinal=counternucleotidefinal+2;
          }
       }
  }
 
 return sentenceAminoOutput;
}


//};


int main()
{
 DNASequence seq;
 
 seq.clear();
 cout<<seq.getLength()<<endl;//prints 0

 seq.clear();
 cout<<seq.getNucleotideSequence()<<endl;//prints "" (nothing)

 seq.appendCodon('G','A','T');
 cout<<seq.getNucleotideSequence()<<endl;//prints GAT
 seq.appendCodon('T','T','X');//nothing is appended, X is invalid
 cout<<seq.getNucleotideSequence()<<endl;//still prints GAT

 seq.appendAminoAcid('Q');//this is Glutamine codon
 cout<<seq.getNucleotideSequence()<<endl;//prints GATCAA

 seq.appendAminoAcid('Z');//Z is invalid; nothing appended
 cout<<seq.getNucleotideSequence()<<endl;//still prints GATCAA

 seq.appendAminoAcidsFromFile("amino_human_iris.txt");//**************************************
 cout<<seq.getLength()<<endl;//prints 15
 cout<<seq.getNucleotideSequence()<<endl;//prints: GATCAATGGGCTTAT

 seq.clear();
 //cout<<seq.getNucleotideSequence()<<endl;//con esto se arregla
 seq.appendCodon('A','T','G');
 //cout<<seq.getNucleotideSequence()<<endl;//ojo agregado ahorita
 seq.appendAminoAcidsFromFile("amino_human_iris.txt");//contains WAZY************************
 seq.appendCodon('T','A','A');
 cout<<seq.getLength()<<endl;//prints out 15
 cout<<seq.getNucleotideSequence()<<endl;//prints out: ATGTGGGCTTATTAA
 cout<<seq.getAminoAcidSequence()<<endl;//prints out: *WAY!



 return 0;
}

