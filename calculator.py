

"""FUMBANI CHIBAKA

BIO 251 Ka/Ks Calculator: FINAL PROJECT
Started:          17th April, 2009
Draft Submitted:  29th April, 2009

"""
#This file contains a program that determines whether the type of selection that a
#given sequence is undergoing.

# ALL METHODS WITH GREEN COMMENTS WERE NOT USED IN THE FINAL CALCULATIONS

# 2 Methods are used to calcuate the Ka/Ks ratio:
# a) The NG Method  (Using Jukes Cantor (JC) models)
# b) The LWL Method ( JC- K2P - K2P )
# c) The MLWL Method (K2P-K2P-K2P)

#REFERENCE: Zhang et al KaKs_Calculator (PAGE. 261)


from support import*
import math


# Extracts codons from a sequence
def getCodons(sequence):
    start= 0;
    array= [];
    for i in range (0, len(sequence)+1,3):
        stop= i
        if i!=0:
            array.append(sequence[start:stop])
        start= i
    return array;


# translateRNA, codonMixer, synoFind, xSite1, xSite2, asssignState, getEvoCode, getSequenceCode were
# used to retrieve NSX assignments for the input sequences

#Key: N= Non-synonymous, S= Synonymous, X= Both Non-Synonymous and Synonymous

#___________________________________________________________________________________________

""" Translates RNA sequence (Output: aminoacid)"""
def translateRNA(sequence):
    aminoSeq=""
    start=0;
    for i in range (3, len(sequence)+1,3):
        stop= i
        aminoSeq+= getAmino (sequence[start:stop])
        start= i
    return aminoSeq

# CodonMixer: Produces all the 18 possible codon combinations
def codonMixer():
    codon= "AUCG";
    array=[];
    for i in range (0,4):
        for j in range (0, 4):
            for k in range (0,4):
                newCodon= codon[i]+ codon[j]+ codon[k];
                array.append(newCodon)
    return array;        
    

""" Synonymous Finder: returns codons that are synonymous to the sequence
of interest"""
def synoFind(codon):
    mixerArray= codonMixer();
    mixerLength = len(mixerArray);

    result=[];
    index=[];
    
    for j in range (0,mixerLength):
        if (translateRNA(codon)== translateRNA(mixerArray[j])):
            if(codon!= mixerArray[j]):
                result.append(mixerArray[j]);
                index.append(j);
    return result; """,index  """ #index used for testing purposes


"""Identifies codon solutions that have both N and S assignement at codon position 1"""
def xSite1(array):
    count=len(array)-1;
    for i in range (0, len(array)-1):
        if ( (array[i])[0:2]== array[i+1][0:2]):
            count-=1;
    return (count);

"""Identifies codon solutions that have both N and S assignement at codon position 3"""
def xSite2(array):
    count=len(array);
    for i in range (0, len(array)-1):
        if ( (array[i])[1:3]== array[i+1][1:3]):
            count-=1;
    return (count);


"""Analyses the codon to find N ans S nucleotide sites"""
def assignState(codon):
    solution=['N','N','N'];
    array= synoFind(codon);
    
    length= len(array);
    
    firstPos = xSite1(array);
    secPos   = xSite2(array);
    
    for i in range (0, length):
        if(codon[0]!= (array[i])[0]):
            solution[0]='S';

        if(codon[2]!= (array[i])[2]):
            solution[2]='S'

    if(firstPos != 0 and firstPos!=3):
        solution[0]='X';

    if (secPos != 0 and secPos!=3):
        solution[2]='X';
        
    return solution;



"""Provides the complete N/S/X code for the sequence"""
def getSequenceCode( sequence):
    solution= [];
    array= getCodons(sequence);
    length= len(array);
    
    for i in range (0, length):
        solution.append (assignState(array[i]))
    return solution;



"""Detects mutations in codons from 2 sequences"""
def getEvoCode(seq1, seq2):
    index="";
    
    array1= (getCodons(seq1));
    array2= (getCodons(seq2));

    for i in range (0, len(array1)):
        if(array1[i]!= array2[i]):
            if (getAmino(array1[i])!= getAmino(array2[i])):
                index+= 'N'
            else:
                index+= 'S'
    return index

#_______________________________________________________________________________________________


#Checks input sequence for error 1 (not divisible by 3)

def clean(sequence):
    
    if (len(sequence)%3 != 0):
        return False
    else:
        return True

#Converts DNA to RNA
def convertDNAtoRNA(sequence):
    result=""
    for i in range(0,len(sequence)):
        if (sequence[i]== 'T'):   
            result += 'U'
        else:
            result += sequence[i]
    return result    

    



# Returns the Synonymity state of a given codon
# Key: N= Non-synonymous, S= Synonymous, X= Both Non-Synonymous and Synonymous
#REFERENCE:  

def getNSX(codon):
    code = {}
    code["F"]  = "NNX"   # Phenylalanine
    code["L"]  = "XNX"   # Leucine
    code["I"]  = "NNX"   # Isoleucine
    code["M"]  = "NNN"   # Methionine
    code["V"]  = "NNS"   # Valine
    code["S1"] = "NNS"   # Serine2        
    code["S2"] = "XNX"   # Serine2
    code["P"]  = "NNS"   # Proline
    code["T"]  = "NNS"   # Threonine
    code["A"]  = "NNS"   # Alanine
    code["Y"]  = "NNX"   # Tyrosine
    code["."]  = "NXX"   # Stop codon: represented by X
    code["H"]  = "NNX"   # Histidine
    code["Q"]  = "NNX"   # Glutamine
    code["N"]  = "NNX"   # Asparagine
    code["K"]  = "NNX"   # Lysine
    code["D"]  = "NNX"   # Aspartic acid
    code["E"]  = "NNX"   # Glutamic acid
    code["C"]  = "NNX"   # Cysteine
    code["W"]  = "NNN"   # Tryptophan
    code["R"]  = "XNX"   # Arginine
    code["G"]  = "NNS"   # Glycine

    return code[ getAmino(codon)]


#returns the degeneracy code for a specified codon
#Key: 0= Non- degenerate, 2= two-fold Degeneracy, 4= FourFold Degeneracy
#REFERENCE: 


def getDeg(codon):
    code = {}
    code["F"]  = "002"   # Phenylalanine
    code["L"]  = "204"   # Leucine
    code["I"]  = "002"   # Isoleucine
    code["M"]  = "000"   # Methionine
    code["V"]  = "004"   # Valine
    code["S1"] = "004"   # Serine2
    code["S2"] = "004"   # Serine2
    code["P"]  = "004"   # Proline
    code["T"]  = "004"   # Threonine
    code["A"]  = "004"   # Alanine
    code["Y"]  = "002"   # Tyrosine
    code["."]  = "022"   # Stop codon: represented by X
    code["H"]  = "002"   # Histidine
    code["Q"]  = "002"   # Glutamine
    code["N"]  = "002"   # Asparagine
    code["K"]  = "002"   # Lysine
    code["D"]  = "002"   # Aspartic acid
    code["E"]  = "002"   # Glutamic acid
    code["C"]  = "002"   # Cysteine
    code["W"]  = "000"   # Tryptophan
    code["R"]  = "204"   # Arginine
    code["G"]  = "004"   # Glycine
    return code[ getAmino(codon)]




#Provides the complete NSX code for a given complete sequence
#Used to analyse both the EVOLVED and ORIGINAL sequences

def getCode(sequence):
    result="";
    index= getCodons(sequence);
    
    for i in range(0,len(index)):
      result+= getNSX(index[i])
    return result;


"""Count the number of N/S/X values decoded from the sequences
def countNSoriginal(inputArray):
    countN=0
    countS=0
    
    for i in range (0, len(inputArray)):
        for j in range (0, len(inputArray[i])):
            if  ((inputArray[i])[j]=='S' or (inputArray[i])[j]=='X' ):
                countS+=1;
            else:
                countN+=1;
    return countN, countS;


def countNSevolution(inputArray):
    countN= 0;
    countS= 0;

    for i in range (0, len(inputArray)):
        if  ((inputArray[i])=='S' or (inputArray[i])=='X' ):
            countS+=1;
        else:
            countN+=1;
    return countN, countS;
    """
#Count the occurance of N and S from the NSX code of the given sequences

def countNSX(code):
    countN = 0;
    countS = 0;

    for i in range (0, len(code)):
        if  (code[i]=='S' or code[i]=='X' ): # Note: X can be both S and N
            countS+=1;
        if  (code[i]=='N' or code[i]=='X' ):
            countN+=1;
    return countN, countS;         


# returns NSX and 024 data for the selected codon
# Used for debugging (comparison between NSX code and 024( Degeneracy code) reveals errors

def getInfo(codon):
    return getNSX(codon), getDeg(codon);

# provides the complete sequence in degeneracy code
def getDegCode(seq1):
    result = "";
    index = getCodons(seq1);
    
    for i in range(0,len(index)):
      result+= getDeg(index[i])
    return result;


# calculates the frequencies of the three different types of degeneracy types (0,2,4)
def count024(code):
    count0 = 0;
    count2 = 0;
    count4 = 0;

    for i in range (0, len(code)):
        if  (code[i]=='0'):
            count0+=1;
        if  (code[i]=='2'):
            count2+=1;
        if  (code[i]=='4'):
            count4+=1;
    return count0, count2, count4;

# calculate average 0,2,4 degeneracy values for two sequences
def average024(seq1,seq2):
    code1= getDegCode(seq1)
    code2= getDegCode(seq2)
    c1,c2,c3= count024(code1)
    d1,d2,d3= count024(code2)
    return (float(c1+d1)/2),(float(c2+d2)/2),(float(c3+d3)/2)


# test if change is Tranversion or Transition
# Key: Pi= TRANSITION   ( A>G or viceversa)
#      Qv= TRANSVERSION ( U>C or viceversa)

def test(p,q):
    if ((p=='A' and q=='G') or (p=='G' and q=='A')):
        return 'Pi'
    if ((p=='U' and q=='C') or (p=='C' and q=='U')):
        return 'Pi' 
    else:
        return 'Qv'


#Calculates the transition/transversion mutation rate ratio

def PiQvRatio(seq1,seq2):
    countPi =0    #transitions
    countQv =0    #transversion
    for i in range (0, len(seq1)):
        check= test(seq1[i],seq2[i])
        if (check== 'Pi'):
            countPi+=1
        else:
            countQv+=1
    return float(countPi)/countQv

# Calculates the Observed transitional (Pi) and transversion (Qv) differences
#REFERENCE:


def getPQvalues (seq1,seq2):
    Pvalues=[0,0,0] #Array elements correspond to site values [L0,L2,L4]
    Qvalues=[0,0,0] #same as above
    degCode1= getDegCode(seq1)

    L0, L2, L4 =average024(seq1,seq2)
    
    for i in range (0, len(seq1)):
        if (seq1[i]!=seq2[i]):
            if ((test(seq1[i],seq2[i]))== 'Pi'):
                code024= degCode1[i]
                if(code024=='0'):
                    calcP= 1/L0
                    Pvalues[0]=calcP
                if(code024=='2'):
                    calcP= 1/L2
                    Pvalues[1]=calcP
                if(code024=='4'):
                    calcP= 1/L4
                    Pvalues[2]=calcP
                    
            if ((test(seq1[i],seq2[i]))== 'Qv'):
                code024= degCode1[i]
                if(code024=='0'):
                    calcQ= 1/L0
                    Qvalues[0]=calcQ
                if(code024=='2'):
                    calcQ= 1/L2
                    Qvalues[1]=calcQ
                if(code024=='4'):
                    calcQ= 1/L4
                    Qvalues[2]=calcQ
        
    return degCode1, L0,L2, L4, Pvalues, Qvalues;
                    


"""MAIN CALCULATIONS"""

# CALCUATION OPTION 1: NG METHOD


# Calculation using the Jukes- Cantor (JC)  Method

def KaKsCalcNG(seq1,seq2):
    seqCode= getCode(seq1);
    evoCode= getEvoCode(seq1,seq2);

    nSeq, sSeq= countNSX(seqCode);
    nEvo, sEvo= countNSX(evoCode);

    Ka= float(sEvo)/sSeq;
    Ks= float(nEvo)/nSeq;

    return Ka,Ks;



# CALCUATION OPTION 3: MLWL METHOD
# Considering the Relative Likehood of nucleotide and Codon Changes


# Calculates Mean and Approximate Error Variance
# REFERENCE: Li, Wu, and Luo P. 152-153

def meanVarianceCalc(degType, seq1, seq2):
    
    degCode1,L0,L2,L4,PArray,QArray = getPQvalues (seq1,seq2)
    
    if(degType==0):
        L=L0
        P=PArray[0]
        Q=QArray[0]
        
    if(degType==2):
        L=L2
        P=PArray[1]
        Q=QArray[1]
        
    if(degType==4):
        L=L4
        P=PArray[2]
        Q=QArray[2]

    ai= 1/(1-(2*P)-Q)
    if(ai<0): #### ??? Didnt know how to deal with negative output (log is taken in next step)
        ai=0.1
    bi= 1/(1-(2*Q))
    if(bi<0): #### ???
        bi=0.1
    ci= (ai-bi)/2
    if(ci<0): #### ???
        ci=0.1      
    di=bi+ci

    #mean of transitional substns per ith site
    Ai=(0.5)*math.log(ai,2.718281828)-(0.25)*math.log(bi,2.718281828);
    #mean of transversional substns per ith site
    Bi=(0.5)*math.log(bi,2.718281828)
    #Approx error variance (transitional substns)
    VAi= ( ((ai*ai)*P +(ci*ci)*Q) - math.pow((ai*P)+(ci*Q),2))/L
    #Approx error variance (transversion substns)
    VBi= (( (bi*bi)*Q )*(1-Q) )/L
    #total number of substitutions per ith site
    Ki=Ai+Bi
    #variance
    VKi = ( ((ai*ai)*P +(di*di)*Q) - math.pow((ai*P)+(ci*Q),2))/L
    
    return Ai,VAi,Bi,VBi,Ki,VKi,L0,L2,L4;


# Calculation of Ka/Ks using Kimuras- two parameter (K2P) model

def KaKsCalcLWL(seq1,seq2):

    A0,VA0,B0,VB0,K0,VK0,L0,L2,L4 = meanVarianceCalc(0, seq1, seq2)
    A2,VA2,B2,VB2,K2,VK2,L0,L2,L4 = meanVarianceCalc(2, seq1, seq2)
    A4,VA4,B4,VB4,K4,VK4,L0,L2,L4 = meanVarianceCalc(4, seq1, seq2)
    
    Ks  = 3*( (L2*A2) + (L4*K4))/(L2 + (3*L4))
    VKs = 9*(L2*L2*VA2 + L4*L4*VK4)/math.pow((L2+ 3*L4),2)
    
    Ka  = 3*(L2*B2 + L0*K0)/(2*L2 + 3*L0)
    VKa = 9*(L2*L2*VB2 + L0*L0*VK0)/math.pow((2*L2+ 3*L0),2)

    return Ka,Ks,VKa,VKs;


# Calculation of Ka/Ks using both JK and K2P models
# Ks- JC
# Ka- K2P
# Corrrecting- K2P

def KaKsCalcMLWL(seq1,seq2):

    k= PiQvRatio(seq1,seq2)
    
    A0,VA0,B0,VB0,K0,VK0,L0,L2,L4 = meanVarianceCalc(0, seq1, seq2)
    A2,VA2,B2,VB2,K2,VK2,L0,L2,L4 = meanVarianceCalc(2, seq1, seq2)
    A4,VA4,B4,VB4,K4,VK4,L0,L2,L4 = meanVarianceCalc(4, seq1, seq2)

    if(k>=2):
        Ka= (L2*B2 + L0*K0)/( ((2*L2)/((k-1)+2))+L0 )
        Ks= (L2*A2 + L4*K4)/(((k-1)*L2/(k-1)+2)+L4)

    else:
        Ka= (L2*B2 + L0*K0)/(((2*L2)/3)+L0 )
        Ks= (L2*A2 + L4*K4)/((L2/3)+L4)

    return Ka,Ks;


    







    
    


        





    
    

                
                
                
    
    
    

         
        
    



    
    

                





            
        
    





