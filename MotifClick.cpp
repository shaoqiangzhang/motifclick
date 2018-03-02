/*===================================================================
Find k-length motifs in the cliques of 2(k-1)-mers (MotifClick)
MotifClick (Version6_BMC2)
Author: Shaoqiang Zhang, March 1 2011
==================================================================*/
#include<map>
#include<string>
#include<iostream>
#include<vector>
#include<fstream>
#include<sstream>
#include<cmath>
#include<ctime>
#include <cstdlib>
//==================================================================
using namespace std;
typedef vector<vector<char> > MatrixChar;
typedef vector<vector<string> > MatrixStr;
typedef vector<vector<int> > MatrixInt;
typedef map<string,int> hashmap;
typedef map<string,vector<int> > vecthash;
typedef map<string,hashmap> hash_of_hash;
//==================================================================
int main(int argc, const char** argv){
    int MotifLength=10;//k-mer
    int MotifNumber=5;//number of top motifs to find
    int densitycutoff=100;//the density cutoff  of the constructed graph
    int sites=1;//simcutoff=simcutoff+sites
    int biDirect=1;//forward and reverse search
    double mergecutoffA=0.5; //mergecutoffA<=mergecutoffB
    double mergecutoffB=0.6;
    double sumSqrDistance=0.3;//Sum of squared distance
    if(argc!=2 && argc!=4 && argc!=6 && argc!=8 && argc!=10 && argc!=12 && argc!=14){
        //cout<<"\nMotifClicker is a program for finding k-length motifs in cliques of 2(k-1)-mers.\n";
        cout<<"\n******\nUSAGE:\n******\n";
        cout<<argv[0]<<" <dataset> [OPTIONS]  > OutputFile\n";
        cout<<"\n<dataset>\tfile containing DNA sequences in FASTA format\n";
        cout<<"\nOPTIONS:\n";
        cout<<"-w\t\tmotif width(default="<<MotifLength<<")\n";
        cout<<"-n\t\tmaximum number of motifs to find(default="<<MotifNumber<<")\n";
		cout<<"-SSD\tupper bound of SSD (sum of squared distance, please select 0<SSD<1, default="<<sumSqrDistance<<")\n";
        cout<<"-b 2\t\tif examine sites on both of DNA strands(default=1 only forward)\n";
        cout<<"-d\t\tupper bound of graph density(default="<<densitycutoff<<")\n";
        cout<<"-s 0\t\tif want more degenerate sites (default=1 if want fewer sites)\n";
	    cout<<"\n*******\nDesigned by Shaoqiang Zhang (zhangshaoqiang@gmail.com), March 2011\n";
        cout<<endl;
        exit(1);
    }
    ifstream fastafile(argv[1]);
    for(int i=2;i<argc-1;i=i+2){
        string strindex(argv[i]);
        if(strindex=="-w"){
            string smotifwidth(argv[i+1]);
            istringstream ismotifwidth(smotifwidth);
            ismotifwidth>>MotifLength;
            if(MotifLength<4){
               cout<<"ERROR: the wrong parameter of '-w', please input an integer(>3)."<<endl; exit(1);
            }
        }else if(strindex=="-n"){
            string smotifnumber(argv[i+1]);
            istringstream ismotifnum(smotifnumber);
            ismotifnum>>MotifNumber;
            if(MotifNumber<1){
                cout<<"ERROR: the wrong parameter of '-n', please input a positive integer."<<endl; exit(1);
            }
		}else if(strindex=="-SSD"){
            string sssd(argv[i+1]);
            istringstream isssd(sssd);
            isssd>>sumSqrDistance;
            if(sumSqrDistance>1 && sumSqrDistance<=0){
                cout<<"ERROR: the wrong parameter of '-SSD', please input 0<SSD<1."<<endl; exit(1);
            }
        }else if(strindex=="-b"){
            string sbidirect(argv[i+1]);
            istringstream isbidirect(sbidirect);
            isbidirect>>biDirect;
            if(biDirect!=1 && biDirect!=2){
                cout<<"ERROR: the wrong parameter of '-b', please input 1 or 2."<<endl; exit(1);
            }
        }else if(strindex=="-s"){
            string ssites(argv[i+1]);
            istringstream issites(ssites);
            issites>>sites;
            if(sites!=1 && sites!=0){
                cout<<"ERROR: the wrong parameter of '-s', please input 0 or 1."<<endl; exit(1);
            }
        }else if(strindex=="-d"){
            string sdensity(argv[i+1]);
            istringstream isdensity(sdensity);
            isdensity>>densitycutoff;
            if(densitycutoff<0){
                cout<<"ERROR: the wrong parameter of '-d', please input an integer(>=0)."<<endl; exit(1);
            }
        }else{
            cout<<"ERROR: the wrong settings of parameters. please type the program name to see the detail!"<<endl; exit(1);
        }
    }
    cout<<"The command line: ";
    for(int j=0;j<argc;j++){
        cout<<argv[j]<<" ";
    }
    cout<<"\n\n";
    //=============================end of input command line===========
    int cnt=0; char ch;
    vector<string> seqName;
    vector<char> seqVector;
    MatrixChar allSeqMat;
    time_t tstart,tend;
    time(&tstart);
    for(string s;getline(fastafile,s);){
        istringstream sin(s);
        sin>>ch;
        if(ch =='>'){
            cnt++;
            seqName.push_back(s);
            if(!seqVector.empty()){
                allSeqMat.push_back(seqVector);
                seqVector.clear();
            }
        }else{
            seqVector.push_back(ch);
            for(char base;sin>>base;){
                seqVector.push_back(base);
            }
        }
    }
    if(cnt==0){
            cout<<"******\nERROR: the input file \"";
            cout<<argv[1]<<"\" isn't FASTA file!!\n******"<<endl;
            exit(1);
    }
    if(!seqVector.empty()){
        allSeqMat.push_back(seqVector);
    }//==================read fasta file into allSeqMat==========================
    int MatSize=allSeqMat.size();
    if(biDirect==2){
        for(int i=0;i<MatSize;++i){
            seqVector.clear();
            for(int j=allSeqMat[i].size()-1;j>=0;j--){
                if(allSeqMat[i][j]=='A' || allSeqMat[i][j]=='a'){
                    seqVector.push_back('T');
                }
                if(allSeqMat[i][j]=='C' || allSeqMat[i][j]=='c'){
                    seqVector.push_back('G');
                }
                if(allSeqMat[i][j]=='G' || allSeqMat[i][j]=='g'){
                    seqVector.push_back('C');
                }
                if(allSeqMat[i][j]=='T' || allSeqMat[i][j]=='t'){
                    seqVector.push_back('A');
                }
            }
            allSeqMat.push_back(seqVector);
        }
    }//=============================push reverse seqs===============================
    int Acnt=0; int Ccnt=0; int Gcnt=0; int Tcnt=0; int allbaseCnt;
    double Abg; double Cbg; double Gbg; double Tbg;
    for(int i=0;i<allSeqMat.size();++i){
        for(int j=0;j<allSeqMat[i].size();++j){
            if(allSeqMat[i][j]=='A' || allSeqMat[i][j]=='a'){
                ++Acnt;
            }
            if(allSeqMat[i][j]=='C' || allSeqMat[i][j]=='c'){
                ++Ccnt;
            }
            if(allSeqMat[i][j]=='G' || allSeqMat[i][j]=='g'){
                ++Gcnt;
            }
            if(allSeqMat[i][j]=='T' || allSeqMat[i][j]=='t'){
                ++Tcnt;
            }
        }
    }

    allbaseCnt=Acnt+Ccnt+Gcnt+Tcnt;
    Abg=static_cast<double>(Acnt)/static_cast<double>(allbaseCnt);
    Cbg=static_cast<double>(Ccnt)/static_cast<double>(allbaseCnt);
    Gbg=static_cast<double>(Gcnt)/static_cast<double>(allbaseCnt);
    Tbg=static_cast<double>(Tcnt)/static_cast<double>(allbaseCnt);
    cout<<"Background distribution of input sequence set:\n";
    cout<<"***********\nA: "<<Abg<<endl;
    cout<<"C: "<<Cbg<<endl;
    cout<<"G: "<<Gbg<<endl;
    cout<<"T: "<<Tbg<<endl;
    cout<<"***********\n"<<endl;
    //double sqrA=pow((Abg-0.5),2)+pow((Tbg-0.5),2)+pow(Cbg,2)+pow(Gbg,2);
    //double sqrB=pow((Cbg-0.5),2)+pow((Gbg-0.5),2)+pow(Abg,2)+pow(Tbg,2);
    //double sumSqrDistance=(sqrA>sqrB)?sqrA:sqrB;//calculate max distance
    //if(sumSqrDistance>0.2){
    //   sumSqrDistance=0.2;
    //}
    //else if(sumSqrDistance<0.16){
    //   sumSqrDistance=0.2;
    //}
    cout<<"SSD cutoff: "<<sumSqrDistance<<"\n"<<endl;//?????????????????
    //=========================calculate background==========
    srand((unsigned)time(0));
    int random_row; int random_column;
    int totalseqnumber=0;
    int allmatchnumber=0;
    int testtimes=allSeqMat.size()/4; //sampling 25% of seqs
    if(testtimes<10){
        testtimes=10;
    }
    for(int lp=0;lp<testtimes;++lp){//test times
        while(1){
            random_row = (rand()%(allSeqMat.size()));
            if(allSeqMat[random_row].size()>=MotifLength){
                random_column=(rand()%(allSeqMat[random_row].size()-MotifLength+1));
                break;
            }
        }
        int currmatchnumber;
        int seqnum=0;
        int maxmatchnumber;
        int totalmatchnumber=0;
        vector<int> sortedMaxMatchVect;sortedMaxMatchVect.clear();
        for(int ind=0;ind<allSeqMat.size();++ind){
            if(allSeqMat[ind].size()<MotifLength ){
                continue;
            }
            maxmatchnumber=0;
            for(int jnd=0;jnd<=allSeqMat[ind].size()-MotifLength;++jnd){
                currmatchnumber=0;
                if(ind==random_row && jnd==random_column){
                        continue;
                }
                for(int len=0;len<MotifLength;++len){
                    if(allSeqMat[ind][jnd+len]==allSeqMat[random_row][random_column+len]){
                         ++currmatchnumber;
                    }
                }
                if(maxmatchnumber<currmatchnumber){
                    maxmatchnumber=currmatchnumber;
                }
            }
            if(sortedMaxMatchVect.empty()){
                sortedMaxMatchVect.push_back(maxmatchnumber);
            }else{
                if(maxmatchnumber>=sortedMaxMatchVect[0]){
                    sortedMaxMatchVect.insert(sortedMaxMatchVect.begin(),1,maxmatchnumber);
                }else if(maxmatchnumber<sortedMaxMatchVect[sortedMaxMatchVect.size()-1]){
                    sortedMaxMatchVect.push_back(maxmatchnumber);
                }else{
                    vector<int>::iterator theIterator = sortedMaxMatchVect.begin();
                    for(int mx=1;mx<sortedMaxMatchVect.size();++mx){
                        ++theIterator;
                        if(maxmatchnumber<sortedMaxMatchVect[mx-1] && maxmatchnumber>=sortedMaxMatchVect[mx]){
                            sortedMaxMatchVect.insert(theIterator,1,maxmatchnumber);
                            break;
                        }
                    }
                }
            }//==========sort by match number==============
        }
        for(int si=0;si<static_cast<int>(static_cast<double>(sortedMaxMatchVect.size()+1)*0.95);++si){
            ++seqnum;
            totalmatchnumber=totalmatchnumber+sortedMaxMatchVect[si];
        }
        totalseqnumber=totalseqnumber+seqnum;
        allmatchnumber=allmatchnumber+totalmatchnumber;
    }
    int averageMatchNumber;
    if(allmatchnumber%totalseqnumber>totalseqnumber/2){//rounding into integers
        averageMatchNumber=allmatchnumber/totalseqnumber+1;
    }else{
        averageMatchNumber=allmatchnumber/totalseqnumber;
    }
    cout<<"Average match number: "<<averageMatchNumber<<"\n"<<endl;
    //=========statistic average match number to decide the cutoff==========================
    int AcntB=0;int CcntB=0;int GcntB=0;int TcntB=0;
    int existA; int existC; int existG; int existT; int existACGT;
    double distA;double distB;double tempDistA; double tempDistB;
    vector<char> SegVectA; vector<char> SegVectB;
    int matchscore;
    int maxMatchNum;int matchNum;
    hash_of_hash nodesHashofHash;
    hashmap subhash;
    subhash.clear(); nodesHashofHash.clear();
    vecthash vhash; vhash.clear();
    vector<int> AMerVect; vector<int> BMerVect;
    string AMerStr; string BMerStr;
    int intervalNumA; int intervalNumB;
    int matchcutoff;
    if(averageMatchNumber==(MotifLength)){
        matchcutoff=MotifLength;
    //}else if(averageMatchNumber==MotifLength){
	//matchcutoff=MotifLength;
    }else{
        matchcutoff=averageMatchNumber+sites;
        //if(MotifLength>10 && static_cast<double>(averageMatchNumber)/static_cast<double>(MotifLength)>=0.8){
	//	matchcutoff=averageMatchNumber+sites-1;
	//}
    }
    //=======================================================
    for(int i=0;i<allSeqMat.size()-1;++i){
        if(allSeqMat[i].size()<MotifLength){
            continue;
        }
        if(allSeqMat[i].size()%(MotifLength-1)==0){
            intervalNumA=(allSeqMat[i].size()/(MotifLength-1))-1;
        }else{
            intervalNumA=(allSeqMat[i].size()/(MotifLength-1));
        }
        for(int j=0;j<intervalNumA;++j){
            SegVectA.clear();
            for(int p=j*(MotifLength-1);p<(j*(MotifLength-1)+2*(MotifLength-1));++p){
                if(p<allSeqMat[i].size()){
                    SegVectA.push_back(allSeqMat[i][p]);
                }
            }
            for(int s=i+1;s<allSeqMat.size();++s){
                if(allSeqMat[s].size()<MotifLength){
                    continue;
                }
                if(allSeqMat[s].size()%(MotifLength-1)==0){
                    intervalNumB=(allSeqMat[s].size()/(MotifLength-1))-1;
                }else{
                    intervalNumB=(allSeqMat[s].size()/(MotifLength-1));
                }
                for(int t=0;t<intervalNumB;++t){
                    SegVectB.clear();
                    for(int q=t*(MotifLength-1);q<(t*(MotifLength-1)+2*(MotifLength-1));++q){
                        if(q<allSeqMat[s].size()){
                            SegVectB.push_back(allSeqMat[s][q]);
                        }
                    }
                    Acnt=Ccnt=Gcnt=Tcnt=0;
                    for(int k=0;k<MotifLength;++k){
                        if(SegVectA[k]=='A' || SegVectA[k]=='a'){
                            Acnt++;
                        }
                        if(SegVectA[k]=='C' || SegVectA[k]=='c'){
                            Ccnt++;
                        }
                        if(SegVectA[k]=='G' || SegVectA[k]=='g'){
                            Gcnt++;
                        }
                        if(SegVectA[k]=='T' || SegVectA[k]=='t'){
                            Tcnt++;
                        }
                    }
                    maxMatchNum=0;
                    for(int a=0;a<=SegVectA.size()-MotifLength;++a){
                        existA=existC=existG=existT=0;
                        if(a>0){
                            if(SegVectA[a-1]=='A' || SegVectA[a-1]=='a'){
                                Acnt--;
                            }
                            if(SegVectA[a+MotifLength-1]=='A' || SegVectA[a+MotifLength-1]=='a'){
                                Acnt++;
                            }
                            if(SegVectA[a-1]=='C' || SegVectA[a-1]=='c'){
                                Ccnt--;
                            }
                            if(SegVectA[a+MotifLength-1]=='C'|| SegVectA[a+MotifLength-1]=='c'){
                                Ccnt++;
                            }
                            if(SegVectA[a-1]=='G' || SegVectA[a-1]=='g'){
                                Gcnt--;
                            }
                            if(SegVectA[a+MotifLength-1]=='G' || SegVectA[a+MotifLength-1]=='g'){
                                Gcnt++;
                            }
                            if(SegVectA[a-1]=='T' || SegVectA[a-1]=='t'){
                                Tcnt--;
                            }
                            if(SegVectA[a+MotifLength-1]=='T' || SegVectA[a+MotifLength-1]=='t'){
                                Tcnt++;
                            }
                        }
                        if(Acnt>0){ existA=1;} if(Ccnt>0){existC=1;} if(Gcnt>0){existG=1;} if(Tcnt>0){existT=1;}
                        existACGT=existA+existC+existG+existT;
                        distA=pow((static_cast<double>(Acnt)/static_cast<double>(MotifLength)-Abg),2)+
                            pow((static_cast<double>(Ccnt)/static_cast<double>(MotifLength)-Cbg),2)+
                            pow((static_cast<double>(Gcnt)/static_cast<double>(MotifLength)-Gbg),2)+
                            pow((static_cast<double>(Tcnt)/static_cast<double>(MotifLength)-Tbg),2);
                        if(distA>=(sumSqrDistance-0.001) || existACGT<3){
                            continue;
                        }else{
                            AcntB=CcntB=GcntB=TcntB=0;
                            for(int k=0;k<MotifLength;++k){
                                if(SegVectB[k]=='A' || SegVectB[k]=='a'){
                                    AcntB++;
                                }
                                if(SegVectB[k]=='C' || SegVectB[k]=='c'){
                                    CcntB++;
                                }
                                if(SegVectB[k]=='G' || SegVectB[k]=='g'){
                                    GcntB++;
                                }
                                if(SegVectB[k]=='T' || SegVectB[k]=='t'){
                                    TcntB++;
                                }
                            }
                            for(int b=0;b<=SegVectB.size()-MotifLength;++b){
                                existA=existC=existG=existT=0;
                                if(b>0){
                                    if(SegVectB[b-1]=='A' || SegVectB[b-1]=='a'){
                                        AcntB--;
                                    }
                                    if(SegVectB[b+MotifLength-1]=='A' || SegVectB[b+MotifLength-1]=='a'){
                                        AcntB++;
                                    }
                                    if(SegVectB[b-1]=='C' || SegVectB[b-1]=='c'){
                                        CcntB--;
                                    }
                                    if(SegVectB[b+MotifLength-1]=='C' || SegVectB[b+MotifLength-1]=='c'){
                                        CcntB++;
                                    }
                                    if(SegVectB[b-1]=='G' || SegVectB[b-1]=='g'){
                                        GcntB--;
                                    }
                                    if(SegVectB[b+MotifLength-1]=='G' || SegVectB[b+MotifLength-1]=='g'){
                                        GcntB++;
                                    }
                                    if(SegVectB[b-1]=='T' || SegVectB[b-1]=='t'){
                                        TcntB--;
                                    }
                                    if(SegVectB[b+MotifLength-1]=='T' || SegVectB[b+MotifLength-1]=='t'){
                                        TcntB++;
                                    }
                                }
                                if(AcntB>0){ existA=1;} if(CcntB>0){existC=1;} if(GcntB>0){existG=1;} if(TcntB>0){existT=1;}
                                existACGT=existA+existC+existG+existT;
                                distB=pow((static_cast<double>(AcntB)/static_cast<double>(MotifLength)-Abg),2)+
                                    pow((static_cast<double>(CcntB)/static_cast<double>(MotifLength)-Cbg),2)+
                                    pow((static_cast<double>(GcntB)/static_cast<double>(MotifLength)-Gbg),2)+
                                    pow((static_cast<double>(TcntB)/static_cast<double>(MotifLength)-Tbg),2);
                                if(distB>=(sumSqrDistance-0.001) || existACGT<3){
                                    continue;
                                }else{
                                    matchNum=0;
                                    for(int k=0;k<MotifLength;++k){
                                        if(SegVectA[a+k] == SegVectB[b+k]){
                                            ++matchNum;
                                        }
                                    }
                                    if(maxMatchNum<matchNum){
                                        maxMatchNum=matchNum;
                                        tempDistA=distA;
                                        tempDistB=distB;
                                    }
                                }
                            }
                        }
                    }
                    matchscore=maxMatchNum;
                    if(matchscore>=matchcutoff){////////////////////////cutoff???
                    //==============the following is similarity graph construction==========
                        stringstream stri; stringstream strj; stringstream strs; stringstream strt;
                        stri<<i; strj<<j; strs<<s; strt<<t;
                        AMerStr=stri.str()+"-"+strj.str();
                        BMerStr=strs.str()+"-"+strt.str();
                        AMerVect.clear(); AMerVect.push_back(i); AMerVect.push_back(j);
                        BMerVect.clear(); BMerVect.push_back(s); BMerVect.push_back(t);
                        if(vhash.count(AMerStr)<1){
                            vhash.insert(vecthash::value_type(AMerStr,AMerVect));
                        }
                        if(vhash.count(BMerStr)<1){
                            vhash.insert(vecthash::value_type(BMerStr,BMerVect));
                        }
                        if(nodesHashofHash.count(AMerStr)<1){
                            subhash.clear();
                            subhash.insert(hashmap::value_type(BMerStr, matchscore));
                            nodesHashofHash.insert(hash_of_hash::value_type(AMerStr,subhash));
                        }else{
                            nodesHashofHash[AMerStr][BMerStr]=matchscore;
                        }
                        if(nodesHashofHash.count(BMerStr)<1){
                            subhash.clear();
                            subhash.insert(hashmap::value_type(AMerStr,matchscore));
                            nodesHashofHash.insert(hash_of_hash::value_type(BMerStr,subhash));
                        }else{
                            nodesHashofHash[BMerStr][AMerStr]=matchscore;
                        }
                    }
                    //gapless local alignment for each pair of 2(k-1)-mers
                }
            }
        }
        //cout<<endl;
    }//==================end of constructing mutipartite similarity graph of 2(k-1)-mers===================

    for(int i=0;i<allSeqMat.size();++i){
        if(allSeqMat[i].size()<MotifLength){
            continue;
        }
        if(allSeqMat[i].size()%(MotifLength-1)==0){
            intervalNumA=(allSeqMat[i].size()/(MotifLength-1))-1;
        }else{
            intervalNumA=(allSeqMat[i].size()/(MotifLength-1));
        }
        for(int j=0;j<intervalNumA-2;++j){
            SegVectA.clear();
            for(int p=j*(MotifLength-1);p<(j*(MotifLength-1)+2*(MotifLength-1));++p){
                if(p<allSeqMat[i].size()){
                    SegVectA.push_back(allSeqMat[i][p]);
                }
            }
            for(int t=j+2;t<intervalNumA;++t){
                SegVectB.clear();
                for(int q=t*(MotifLength-1);q<(t*(MotifLength-1)+2*(MotifLength-1));++q){
                    if(q<allSeqMat[i].size()){
                        SegVectB.push_back(allSeqMat[i][q]);
                    }
                }
                Acnt=Ccnt=Gcnt=Tcnt=0;
                for(int k=0;k<MotifLength;++k){
                    if(SegVectA[k]=='A' || SegVectA[k]=='a'){
                        Acnt++;
                    }
                    if(SegVectA[k]=='C' || SegVectA[k]=='c'){
                        Ccnt++;
                    }
                    if(SegVectA[k]=='G' || SegVectA[k]=='g'){
                        Gcnt++;
                    }
                    if(SegVectA[k]=='T' || SegVectA[k]=='t'){
                        Tcnt++;
                    }
                }
                maxMatchNum=0;
                for(int a=0;a<=SegVectA.size()-MotifLength;++a){
                    existA=existC=existG=existT=0;
                    if(a>0){
                        if(SegVectA[a-1]=='A' || SegVectA[a-1]=='a'){
                            Acnt--;
                        }
                        if(SegVectA[a+MotifLength-1]=='A' || SegVectA[a+MotifLength-1]=='a'){
                            Acnt++;
                        }
                        if(SegVectA[a-1]=='C' || SegVectA[a-1]=='c'){
                            Ccnt--;
                        }
                        if(SegVectA[a+MotifLength-1]=='C'|| SegVectA[a+MotifLength-1]=='c'){
                            Ccnt++;
                        }
                        if(SegVectA[a-1]=='G' || SegVectA[a-1]=='g'){
                            Gcnt--;
                        }
                        if(SegVectA[a+MotifLength-1]=='G' || SegVectA[a+MotifLength-1]=='g'){
                            Gcnt++;
                        }
                        if(SegVectA[a-1]=='T' || SegVectA[a-1]=='t'){
                            Tcnt--;
                        }
                        if(SegVectA[a+MotifLength-1]=='T' || SegVectA[a+MotifLength-1]=='t'){
                            Tcnt++;
                        }
                    }
                    if(Acnt>0){ existA=1;} if(Ccnt>0){existC=1;} if(Gcnt>0){existG=1;} if(Tcnt>0){existT=1;}
                    existACGT=existA+existC+existG+existT;
                    distA=pow((static_cast<double>(Acnt)/static_cast<double>(MotifLength)-Abg),2)+
                        pow((static_cast<double>(Ccnt)/static_cast<double>(MotifLength)-Cbg),2)+
                        pow((static_cast<double>(Gcnt)/static_cast<double>(MotifLength)-Gbg),2)+
                        pow((static_cast<double>(Tcnt)/static_cast<double>(MotifLength)-Tbg),2);
                    if(distA>=(sumSqrDistance-0.001) || existACGT<3){
                        continue;
                    }else{
                        AcntB=CcntB=GcntB=TcntB=0;
                        for(int k=0;k<MotifLength;++k){
                            if(SegVectB[k]=='A' || SegVectB[k]=='a'){
                                AcntB++;
                            }
                            if(SegVectB[k]=='C' || SegVectB[k]=='c'){
                                CcntB++;
                            }
                            if(SegVectB[k]=='G' || SegVectB[k]=='g'){
                                GcntB++;
                            }
                            if(SegVectB[k]=='T' || SegVectB[k]=='t'){
                                TcntB++;
                            }
                        }
                        for(int b=0;b<=SegVectB.size()-MotifLength;++b){
                            existA=existC=existG=existT=0;
                            if(b>0){
                                if(SegVectB[b-1]=='A' || SegVectB[b-1]=='a'){
                                    AcntB--;
                                }
                                if(SegVectB[b+MotifLength-1]=='A' || SegVectB[b+MotifLength-1]=='a'){
                                    AcntB++;
                                }
                                if(SegVectB[b-1]=='C' || SegVectB[b-1]=='c'){
                                    CcntB--;
                                }
                                if(SegVectB[b+MotifLength-1]=='C' || SegVectB[b+MotifLength-1]=='c'){
                                    CcntB++;
                                }
                                if(SegVectB[b-1]=='G' || SegVectB[b-1]=='g'){
                                    GcntB--;
                                }
                                if(SegVectB[b+MotifLength-1]=='G' || SegVectB[b+MotifLength-1]=='g'){
                                    GcntB++;
                                }
                                if(SegVectB[b-1]=='T' || SegVectB[b-1]=='t'){
                                    TcntB--;
                                }
                                if(SegVectB[b+MotifLength-1]=='T' || SegVectB[b+MotifLength-1]=='t'){
                                    TcntB++;
                                }
                            }
                            if(AcntB>0){ existA=1;} if(CcntB>0){existC=1;} if(GcntB>0){existG=1;} if(TcntB>0){existT=1;}
                            existACGT=existA+existC+existG+existT;
                            distB=pow((static_cast<double>(AcntB)/static_cast<double>(MotifLength)-Abg),2)+
                                pow((static_cast<double>(CcntB)/static_cast<double>(MotifLength)-Cbg),2)+
                                pow((static_cast<double>(GcntB)/static_cast<double>(MotifLength)-Gbg),2)+
                                pow((static_cast<double>(TcntB)/static_cast<double>(MotifLength)-Tbg),2);
                            if(distB>=(sumSqrDistance-0.001) || existACGT<3){
                                continue;
                            }else{
                                matchNum=0;
                                for(int k=0;k<MotifLength;++k){
                                    if(SegVectA[a+k] == SegVectB[b+k]){
                                        ++matchNum;
                                    }
                                }
                                if(maxMatchNum<matchNum){
                                    maxMatchNum=matchNum;
                                    tempDistA=distA;
                                    tempDistB=distB;
                                }
                            }
                        }
                    }
                }
                matchscore=maxMatchNum;
                if(matchscore>=matchcutoff){////////////////////////cutoff???
                //==============the following is similarity graph construction==========
                    stringstream stri; stringstream strj; stringstream strt;
                    stri<<i; strj<<j; strt<<t;
                    AMerStr=stri.str()+"-"+strj.str();
                    BMerStr=stri.str()+"-"+strt.str();
                    AMerVect.clear(); AMerVect.push_back(i); AMerVect.push_back(j);
                    BMerVect.clear(); BMerVect.push_back(i); BMerVect.push_back(t);
                    if(vhash.count(AMerStr)<1){
                        vhash.insert(vecthash::value_type(AMerStr,AMerVect));
                    }
                    if(vhash.count(BMerStr)<1){
                        vhash.insert(vecthash::value_type(BMerStr,BMerVect));
                    }
                    if(nodesHashofHash.count(AMerStr)<1){
                        subhash.clear();
                        subhash.insert(hashmap::value_type(BMerStr, matchscore));
                        nodesHashofHash.insert(hash_of_hash::value_type(AMerStr,subhash));
                    }else{
                        nodesHashofHash[AMerStr][BMerStr]=matchscore;
                    }
                    if(nodesHashofHash.count(BMerStr)<1){
                        subhash.clear();
                        subhash.insert(hashmap::value_type(AMerStr,matchscore));
                        nodesHashofHash.insert(hash_of_hash::value_type(BMerStr,subhash));
                    }else{
                        nodesHashofHash[BMerStr][AMerStr]=matchscore;
                    }
                }
                    //gapless local alignment for each pair of 2(k-1)-mers
            }
        }
        //cout<<endl;
    }//=======construct the edges in  the same sequence============================

    int graphcount=0;
    while(1){
        int edgeNum=0;
        for (hash_of_hash::const_iterator hh=nodesHashofHash.begin(); hh!=nodesHashofHash.end(); ++hh){
                edgeNum=edgeNum+(hh->second.size());
        }
        graphcount++;
        int graphdensity=edgeNum/(2*nodesHashofHash.size());
        cout<<"Graph#"<<graphcount<<" density: "<<graphdensity<<" (Match cutoff: "<<matchcutoff<<")\n"<<endl;
        if(graphdensity>densitycutoff){//graph density cutoff???????????????????????
            if(matchcutoff>(MotifLength-1)){
                cout<<"WARNING: We cannot get a graph with density <"<<densitycutoff<<"\n"<<endl;
                break;//exit(1);
            }else{
                hash_of_hash tempNodeHashofHash;
                for (hash_of_hash::const_iterator hh=nodesHashofHash.begin(); hh!=nodesHashofHash.end(); ++hh){
                    hashmap tempNeighborNodeMap=hh->second;
                    for(hashmap::const_iterator h=hh->second.begin();h!=hh->second.end();++h){
                        if(h->second ==matchcutoff){
                            tempNeighborNodeMap.erase(h->first);
                        }
                    }
                    if(!tempNeighborNodeMap.empty()){
                        tempNodeHashofHash.insert(hash_of_hash::value_type(hh->first,tempNeighborNodeMap));
                    }
                }
                nodesHashofHash.clear();
                nodesHashofHash=tempNodeHashofHash;
                tempNodeHashofHash.clear();
                matchcutoff=matchcutoff+1;
            }
        }else{
            break;
        }
    }
    //========calculate the density of the graph and update the graph if it has high density=======
    string currentNode;
    hashmap neighborNodeHash;
    int localNeighborNum;
    string currentNeighborNode;
    hashmap currentNeighborHash;
    int minDegree;
    string minDegreeNode;
    int removeNum;
    int currentEdgesWeightSum;
    int minWeight;
    MatrixStr matStr;//store all cliques, each row is a clique
    //--------------------------------------------------------------
    for (hash_of_hash::const_iterator hh=nodesHashofHash.begin(); hh!=nodesHashofHash.end(); ++hh){
        currentNode=hh->first;
        neighborNodeHash=hh->second;
        if(neighborNodeHash.size()>=2 ){ // && neighborNodeHash.size()>=graphdensity){//??????????????????????????????
            vector<string> currentClique;
            currentClique.push_back(currentNode);
            //cout<<currentNode;
            while(1){
                minDegree=neighborNodeHash.size();
                minWeight=minDegree;
                for(hashmap::const_iterator h=neighborNodeHash.begin();h!=neighborNodeHash.end();++h){
                    //cout<<"("<<h->first<<" -> "<<h->second<<")";
                    localNeighborNum=1;
                    currentNeighborNode=h->first;
                    currentEdgesWeightSum=h->second;
                    hash_of_hash::iterator chh=nodesHashofHash.find(currentNeighborNode);
                    //cout<<chh->first;
                    currentNeighborHash=chh->second;
                    for(hashmap::const_iterator ch=currentNeighborHash.begin();ch!=currentNeighborHash.end();++ch){
                        if(neighborNodeHash.count(ch->first)>0){
                            localNeighborNum++;
                            currentEdgesWeightSum=+ch->second;
                        }
                    }
                    //cout<<"degree: "<<localNeighborNum<<" ) ";
                    if(localNeighborNum<minDegree){
                        minDegree=localNeighborNum;
                        minDegreeNode=currentNeighborNode;
                        minWeight=currentEdgesWeightSum;
                    }else if(localNeighborNum==minDegree){
                        if(currentEdgesWeightSum<minWeight){
                            minDegree=localNeighborNum;
                            minDegreeNode=currentNeighborNode;
                            minWeight=currentEdgesWeightSum;
                        }
                    }

                }
                if(minDegree == neighborNodeHash.size()){
                    break;
                }else{
                    removeNum=neighborNodeHash.erase(minDegreeNode);
                }
            }
            for(hashmap::const_iterator ih=neighborNodeHash.begin();ih!=neighborNodeHash.end();++ih){
                currentClique.push_back(ih->first);
                //cout<<" "<<ih->first;
            }
            matStr.push_back(currentClique);
            //cout<<endl;
        }
    }//===========find the max clique associated with each node===============
    vector<int> cliquesWeightVect;
    hashmap cneighborhash;
    for(int ci=0;ci<matStr.size();++ci){
        int cliqWeightSum=0;
        for(int cj=0;cj<matStr[ci].size()-1;++cj){
            //cout<<matStr[ci][cj]<<" ";
            hash_of_hash::iterator cnode=nodesHashofHash.find(matStr[ci][cj]);
            cneighborhash=cnode->second;
            for(int ck=cj+1;ck<matStr[ci].size();++ck){
                hashmap::iterator cn=cneighborhash.find(matStr[ci][ck]);
                cliqWeightSum=cliqWeightSum+(cn->second);
            }
        }
        //cout<<matStr[ci][matStr[ci].size()-1]<<" ";
        cliquesWeightVect.push_back(cliqWeightSum);
        //cout<<cliqWeightSum<<endl;
    }
    //========calculate the sum of weights for each clique===================
    int maxWeight;
    int maxWeightIndex;
    vector<int> MergedIndex_of_matStr;
    vector<string> cliqueVect;
    vector<string> quasicliqueVect;
    vector<string> unmatchVect;
    MatrixStr tempCliquesMatStr;
    vector<int> tempCliquesWeightVect;
    for(int lp=1;lp<=MotifNumber;++lp){
        //===============find top number of motifs====
        maxWeight=0;
        maxWeightIndex;
        if(cliquesWeightVect.size()==0){
            break;
        }else{
            for(int w=0;w<cliquesWeightVect.size();++w){
                if(maxWeight<cliquesWeightVect[w]){
                    maxWeight=cliquesWeightVect[w];
                    maxWeightIndex=w;
                }
            }
        }
        cout<<"\n\n*******\nMOTIF#"<<lp<<endl;
        cout<<"The corresponding clique (Sum of weights="<<maxWeight<<"): ";
        for(int z=0;z<matStr[maxWeightIndex].size();++z){
            cout<<" "<<matStr[maxWeightIndex][z];
        }
        cout<<endl;
        //======find the best clique by sorting weights========================
        MergedIndex_of_matStr.clear();
        //MergedIndex_of_matStr.push_back(maxWeightIndex);
        cliqueVect=matStr[maxWeightIndex];
        quasicliqueVect=matStr[maxWeightIndex];
        matStr.erase(matStr.begin()+maxWeightIndex);
        cliquesWeightVect.erase(cliquesWeightVect.begin()+maxWeightIndex);
        for(int mi=0;mi<matStr.size(); ++mi){
            //if(static_cast<double>(matStr[mi].size())/static_cast<double>(matStr[maxWeightIndex].size())>=0.25){
                int matchCount=0;
                unmatchVect.clear();
                for(int mj=0;mj<matStr[mi].size();++mj){
                    bool judgeExistInMaxWeight=0;
                    for(int mm=0;mm<cliqueVect.size();++mm){
                        if(cliqueVect[mm]==matStr[mi][mj]){
                            matchCount++;
                            break;
                        }
                    }
                    for(int mq=0;mq<quasicliqueVect.size();++mq){
                        if(quasicliqueVect[mq]==matStr[mi][mj]){
                            judgeExistInMaxWeight=1;
                            break;
                        }
                    }
                    if(!judgeExistInMaxWeight){
                        unmatchVect.push_back(matStr[mi][mj]);
                    }
                }
                if(static_cast<double>(matchCount)/static_cast<double>(matStr[maxWeightIndex].size())>=mergecutoffA
                  || static_cast<double>(matchCount)/static_cast<double>(matStr[mi].size())>=mergecutoffB){//merge cliques' cutoff???????????????????
                    MergedIndex_of_matStr.push_back(mi);
                    for(int mp=0;mp<unmatchVect.size();++mp){
                        quasicliqueVect.push_back(unmatchVect[mp]);
                    }
                }
            //}
        }//===========merge the best clique with other cliques into a quasi-clique=======
        tempCliquesWeightVect.clear();
        tempCliquesMatStr.clear();
        for(int ms=0;ms<matStr.size();++ms){
            bool bl=0;
            for(int mn=0;mn<MergedIndex_of_matStr.size();++mn){
                if(ms == MergedIndex_of_matStr[mn]){
                    bl=1;
                    break;
                }
            }
            if(!bl){
                tempCliquesMatStr.push_back(matStr[ms]);
                tempCliquesWeightVect.push_back(cliquesWeightVect[ms]);
            }
        }
        matStr=tempCliquesMatStr;
        cliquesWeightVect=tempCliquesWeightVect;
        //==========update the vector of cliques (matStr)=======================
        vector<int> nodeIndexVect;
        vector<char> segVect;
        MatrixChar topCliqMat;
        MatrixInt topCliqStartPositionMat;
        topCliqMat.clear();
        topCliqStartPositionMat.clear();
        cout<<"Expand it into a quasi-clique: ";
        for(int y=0;y<quasicliqueVect.size();++y){
            cout<<" "<<quasicliqueVect[y];
            vecthash::iterator vh=vhash.find(quasicliqueVect[y]);
            nodeIndexVect=vh->second;
            topCliqStartPositionMat.push_back(nodeIndexVect);
            int ii=nodeIndexVect[0]; int jj=nodeIndexVect[1];
            //cout<<seqName[ii]<<"\n"<<jj*(MotifLength-1)+1<<"\t"<<jj*(MotifLength-1)+2*(MotifLength-1)<<"\t";
            segVect.clear();
            for(int p=jj*(MotifLength-1);p<(jj*(MotifLength-1)+2*(MotifLength-1));++p){
                if(p<allSeqMat[ii].size()){
                    segVect.push_back(allSeqMat[ii][p]);
                    //cout<<allSeqMat[ii][p];
                }
            }
            topCliqMat.push_back(segVect);
            //cout<<endl;
        }
        cout<<"\nDo local gapless alignments and ouput the final motif as follows:\n\n";
        //======store the quasiclique in Matrix topCliqMat===============
        MatrixInt BestAlignVect;
        MatrixInt CurrentAlignVect;
        BestAlignVect.clear();
        int maxtotalcount=0;
        for(int f=0;f<topCliqMat.size();++f){
            for(int k=0;k<=topCliqMat[f].size()-MotifLength;++k){
                CurrentAlignVect.clear();
                vector<int> currentIndex;
                currentIndex.clear();
                currentIndex.push_back(f);currentIndex.push_back(k);
                CurrentAlignVect.push_back(currentIndex);
                int totalcount=0;
                for(int s=0;s<topCliqMat.size() && s!=f;++s){
                    int maxseqcount=0;
                    int maxt=0;
                    for(int t=0;t<=topCliqMat[s].size()-MotifLength;++t){
                        int seqcount=0;
                        for(int r=0;r<MotifLength;++r){
                            if(topCliqMat[f][k+r] ==topCliqMat[s][t+r]){
                                seqcount++;
                            }
                        }
                        if(maxseqcount<seqcount){
                            maxseqcount=seqcount;
                            maxt=t;
                        }
                        if(maxseqcount == MotifLength){
                            break;
                        }
                    }
                    if(maxseqcount>=(averageMatchNumber)){//averageMatchNumber or matchcutoff or matchcutoff-1 ??????
                        totalcount=totalcount+maxseqcount;
                        currentIndex.clear();
                        currentIndex.push_back(s);currentIndex.push_back(maxt);
                        CurrentAlignVect.push_back(currentIndex);
                    }
                }
                if(BestAlignVect.size()<CurrentAlignVect.size()){
                    maxtotalcount=totalcount;
                    BestAlignVect=CurrentAlignVect;
                }else if (BestAlignVect.size()==CurrentAlignVect.size()){
                    if(maxtotalcount<totalcount){
                        maxtotalcount=totalcount;
                        BestAlignVect=CurrentAlignVect;
                    }
                }
            }
        }
        for(int f=0;f<BestAlignVect.size();++f){
            if(biDirect==1){
                cout<<seqName[topCliqStartPositionMat[BestAlignVect[f][0]][0]]<<"\tStart_position: ";
                cout<<topCliqStartPositionMat[BestAlignVect[f][0]][1]*(MotifLength-1)+1+BestAlignVect[f][1]<<endl;
            }
            if(biDirect==2){
                if(topCliqStartPositionMat[BestAlignVect[f][0]][0]<MatSize){
                    cout<<seqName[topCliqStartPositionMat[BestAlignVect[f][0]][0]]<<"\tForward\tStart_position: ";
                    cout<<topCliqStartPositionMat[BestAlignVect[f][0]][1]*(MotifLength-1)+1+BestAlignVect[f][1]<<endl;
                }else{
                    cout<<seqName[topCliqStartPositionMat[BestAlignVect[f][0]][0]-MatSize]<<"\tReverse\tStart_position: ";
                    cout<<allSeqMat[topCliqStartPositionMat[BestAlignVect[f][0]][0]].size()-
                    (topCliqStartPositionMat[BestAlignVect[f][0]][1]*(MotifLength-1)+1+BestAlignVect[f][1])+1<<endl;
                }
            }
            for(int t=BestAlignVect[f][1];t<BestAlignVect[f][1]+MotifLength;++t){
                cout<<topCliqMat[BestAlignVect[f][0]][t];
            }
            cout<<endl;
        }
        //======do alignments for the top clique==========================
    }//=====finish the top number of motifs prediction
    time(&tend);
    double dif = difftime (tend,tstart);
    cout<<"\n\nTotal running time: "<<dif<<" seconds\n"<<endl;
}

