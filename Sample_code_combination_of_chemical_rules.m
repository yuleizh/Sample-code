clear;clc;
%Setall=xlsread('combinedflags_first_2','Sheet1','AC4:AH4');
filename='singleflags_v2';
Total=xlsread(filename,'clinical','A1:AU138');
PSR=xlsread(filename,'PSR','A1:AU138');
ACSINS=xlsread(filename,'ACSINS','A1:AU138');
CSI=xlsread(filename,'CSI','A1:AU138');
ELISA=xlsread(filename,'ELISA','A1:AU138');
BVP=xlsread(filename,'BVP','A1:AU138');
test=xlsread(filename,'test_adi','A1:AU500');
train=xlsread(filename,'train_adi','A1:AU500');
Set=[2 11 15 31 37];
Total1=(sum(Total(:,Set),2)>=2)*1;
PSR1=sum(PSR(:,Set),2)>=2*1;
ACSINS1=sum(ACSINS(:,Set),2)>=2*1;
CSI1=sum(CSI(:,Set),2)>=2*1;
ELISA1=sum(ELISA(:,Set),2)>=2*1;
BVP1=sum(BVP(:,Set),2)>=2*1;
test1=sum(test(:,Set),2)>=2*1;
train1=(sum(train(:,Set),2)>=2)*1;
Total2=find(Total1==0);PSR2=find(PSR1==0);ACSINS2=find(ACSINS1==0);CSI2=find(CSI1==0);
ELISA2=find(ELISA1==0);BVP2=find(BVP1==0);train2=find(train1==0);
hc=sum(Total2<=97);lc=sum(Total2>97);ha=sum(train2<=281);la=sum(train2>281);
gp=sum(PSR2<=109);bp=sum(PSR2>109);ga=sum(ACSINS2<=109);ba=sum(ACSINS2>109);
gc=sum(CSI2<=109);bc=sum(CSI2>109);ge=sum(ELISA2<=107);be=sum(ELISA2>107);gb=sum(BVP2<=102);bb=sum(BVP2>102);
Finalflag=[];
high=10;low=20;
Regionx={'CDR','VH','Fv','H1','H2','H3','L1','L2','L3','VL'};
Finalflag=[];
for region=1:10
    clear opt opt_0 opt_1 opt_2 Flag_AA1 Sign1 Sign_v Res_1 Row
Region=Regionx{region};
warning('off','all');
CAB=xlsread('Clinical_AA',Region,'A1:T200');
%HAB=xlsread('human antibody_sequence_AA comp_newlimit_v1',Region,'A2:T20013');
CP=xlsread('Clinical_AA_PSR',Region,'A1:T300');
CA=xlsread('Clinical_AA_ACSINS',Region,'A1:T300');
CC=xlsread('Clinical_AA_CSI',Region,'A1:T300');
CE=xlsread('Clinical_AA_ELISA',Region,'A1:T300');
CB=xlsread('Clinical_AA_BVP',Region,'A1:T300');
AAB=xlsread('Human_AA_100sim',Region,'A1:T2000');

adihs=[];adils=[];clis_ph=[];clis_pl=[];clis_ch=[];clis_cl=[];clis_eh=[];clis_el=[];clis_bh=[];clis_bl=[];
clis_ah=[];clis_al=[];clihs=[];clils=[];

adi=xlsread('2ndset5_combine_singleflags_v2_new','adimab','A2:CV100');
clis_p=xlsread('2ndset5_combine_singleflags_v2_new','PSR','A2:CV20');
clis_a=xlsread('2ndset5_combine_singleflags_v2_new','ACSINS','A2:CV20');
clis_c=xlsread('2ndset5_combine_singleflags_v2_new','CSI','A2:CV20');
clis_e=xlsread('2ndset5_combine_singleflags_v2_new','ELISA','A2:CV20');
clis_b=xlsread('2ndset5_combine_singleflags_v2_new','BVP','A2:CV20');
clim=xlsread('2ndset5_combine_singleflags_v2_new','clinical','A2:CV100');
ngp=sum(clis_p(:,1)<=gp);nbp=sum(clis_p(:,1)>gp);
nga=sum(clis_a(:,1)<=ga);nba=sum(clis_a(:,1)>ga);
ngc=sum(clis_c(:,1)<=gc);nbc=sum(clis_c(:,1)>gc);
nge=sum(clis_e(:,1)<=ge);nbe=sum(clis_e(:,1)>ge);
ngb=sum(clis_b(:,1)<=gb);nbb=sum(clis_b(:,1)>gb);
nha=sum(adi(:,1)<=ha);nla=sum(adi(:,1)>ha);
nhc=sum(clim(:,1)<=hc);nlc=sum(clim(:,1)>hc);


%hv=1033;lv=59;hr=104;lr=6;highP=0.15;
A=CAB(Total2,1);C=CAB(Total2,2);D=CAB(Total2,3);E=CAB(Total2,4);F=CAB(Total2,5);G=CAB(Total2,6);H=CAB(Total2,7);
I=CAB(Total2,8);K=CAB(Total2,9);L=CAB(Total2,10);M=CAB(Total2,11);N=CAB(Total2,12);P=CAB(Total2,13);Q=CAB(Total2,14);
R=CAB(Total2,15);S=CAB(Total2,16);T=CAB(Total2,17);V=CAB(Total2,18);W=CAB(Total2,19);Y=CAB(Total2,20);
netcharge=R+K+0.1*H-D-E;
max_nc=floor(max(netcharge));
%% Human Antibody information PSR
A1=CP(PSR2,1);C1=CP(PSR2,2);D1=CP(PSR2,3);E1=CP(PSR2,4);F1=CP(PSR2,5);G1=CP(PSR2,6);H1=CP(PSR2,7);
I1=CP(PSR2,8);K1=CP(PSR2,9);L1=CP(PSR2,10);M1=CP(PSR2,11);N1=CP(PSR2,12);P1=CP(PSR2,13);Q1=CP(PSR2,14);
R1=CP(PSR2,15);S1=CP(PSR2,16);T1=CP(PSR2,17);V1=CP(PSR2,18);W1=CP(PSR2,19);Y1=CP(PSR2,20);
Pnetcharge=R1+K1+0.1*H1-D1-E1;
%% Human Antibody information ACSINS
A2=CA(ACSINS2,1);C2=CA(ACSINS2,2);D2=CA(ACSINS2,3);E2=CA(ACSINS2,4);F2=CA(ACSINS2,5);G2=CA(ACSINS2,6);H2=CA(ACSINS2,7);
I2=CA(ACSINS2,8);K2=CA(ACSINS2,9);L2=CA(ACSINS2,10);M2=CA(ACSINS2,11);N2=CA(ACSINS2,12);P2=CA(ACSINS2,13);Q2=CA(ACSINS2,14);
R2=CA(ACSINS2,15);S2=CA(ACSINS2,16);T2=CA(ACSINS2,17);V2=CA(ACSINS2,18);W2=CA(ACSINS2,19);Y2=CA(ACSINS2,20);
Anetcharge=R2+K2+0.1*H2-D2-E2;
%% Human Antibody information CSI
A3=CC(CSI2,1);D3=CC(CSI2,3);E3=CC(CSI2,4);F3=CC(CSI2,5);G3=CC(CSI2,6);H3=CC(CSI2,7);
I3=CC(CSI2,8);K3=CC(CSI2,9);L3=CC(CSI2,10);M3=CC(CSI2,11);N3=CC(CSI2,12);P3=CC(CSI2,13);Q3=CC(CSI2,14);
R3=CC(CSI2,15);S3=CC(CSI2,16);T3=CC(CSI2,17);V3=CC(CSI2,18);W3=CC(CSI2,19);Y3=CC(CSI2,20);
Cnetcharge=R3+K3+0.1*H3-D3-E3;
%% Human Antibody information ELISA
A4=CE(ELISA2,1);C4=CE(ELISA2,2);D4=CE(ELISA2,3);E4=CE(ELISA2,4);F4=CE(ELISA2,5);G4=CE(ELISA2,6);H4=CE(ELISA2,7);
I4=CE(ELISA2,8);K4=CE(ELISA2,9);L4=CE(ELISA2,10);M4=CE(ELISA2,11);N4=CE(ELISA2,12);P4=CE(ELISA2,13);Q4=CE(ELISA2,14);
R4=CE(ELISA2,15);S4=CE(ELISA2,16);T4=CE(ELISA2,17);V4=CE(ELISA2,18);W4=CE(ELISA2,19);Y4=CE(ELISA2,20);
Enetcharge=R4+K4+0.1*H4-D4-E4;
%% Human Antibody information ELISA
A5=CB(BVP2,1);C5=CB(BVP2,2);D5=CB(BVP2,3);E5=CB(BVP2,4);F5=CB(BVP2,5);G5=CB(BVP2,6);H5=CB(BVP2,7);
I5=CB(BVP2,8);K5=CB(BVP2,9);L5=CB(BVP2,10);M5=CB(BVP2,11);N5=CB(BVP2,12);P5=CB(BVP2,13);Q5=CB(BVP2,14);
R5=CB(BVP2,15);S5=CB(BVP2,16);T5=CB(BVP2,17);V5=CB(BVP2,18);W5=CB(BVP2,19);Y5=CB(BVP2,20);
Bnetcharge=R5+K5+0.1*H5-D5-E5;
%% Human B cell from Adimab
ABA=AAB(train2,1);ABC=AAB(train2,2);ABD=AAB(train2,3);ABE=AAB(train2,4);ABF=AAB(train2,5);ABG=AAB(train2,6);ABH=AAB(train2,7);
ABI=AAB(train2,8);ABK=AAB(train2,9);ABL=AAB(train2,10);ABM=AAB(train2,11);ABN=AAB(train2,12);ABP=AAB(train2,13);ABQ=AAB(train2,14);
ABR=AAB(train2,15);ABS=AAB(train2,16);ABT=AAB(train2,17);ABV=AAB(train2,18);ABW=AAB(train2,19);ABY=AAB(train2,20);
ABnetcharge=ABR+ABK+0.1*ABH-ABD-ABE;
%% combination
% hydrophybic residue GAVLIMFWPYRKH
Hydro=[G,A,V,L,I,M,F,W,P,Y,R,K,H];
PHydro=[G1,A1,V1,L1,I1,M1,F1,W1,P1,Y1,R1,K1,H1];
AHydro=[G2,A2,V2,L2,I2,M2,F2,W2,P2,Y2,R2,K2,H2];
CHydro=[G3,A3,V3,L3,I3,M3,F3,W3,P3,Y3,R3,K3,H3];
EHydro=[G4,A4,V4,L4,I4,M4,F4,W4,P4,Y4,R4,K4,H4];
BHydro=[G5,A5,V5,L5,I5,M5,F5,W5,P5,Y5,R5,K5,H5];
AdHydro=[ABG,ABA,ABV,ABL,ABI,ABM,ABF,ABW,ABP,ABY,ABR,ABK,ABH];
Hydro_num=[1,2,3,4,5,6,7,8,9,10,11,12,13];
for i=1:13
    Hydro1{i}=combntns(1:13,i);
    [r(i),~]=size(Hydro1{i});
end    
for k=1:13
    for i=1:1:r(k)
        Hydro_1{k}(:,i)=sum(Hydro(:,Hydro1{k}(i,:)),2);
        PHydro_1{k}(:,i)=sum(PHydro(:,Hydro1{k}(i,:)),2);
        AHydro_1{k}(:,i)=sum(AHydro(:,Hydro1{k}(i,:)),2);
        CHydro_1{k}(:,i)=sum(CHydro(:,Hydro1{k}(i,:)),2);
        EHydro_1{k}(:,i)=sum(EHydro(:,Hydro1{k}(i,:)),2);
        BHydro_1{k}(:,i)=sum(BHydro(:,Hydro1{k}(i,:)),2);
        AdHydro_1{k}(:,i)=sum(AdHydro(:,Hydro1{k}(i,:)),2);
    end
    [~,r_1(k)]=size(Hydro_1{k});
    for j=1:1:r_1(k)
        maxbind_hy{k}(j)=max(Hydro_1{k}(:,j));
        minbind_hy{k}(j)=min(Hydro_1{k}(:,j));
        hy{k}(j)=(maxbind_hy{k}(j)-minbind_hy{k}(j))+1;
    end        
end    
for x=1:1:13
    q1=1;
    opt_0=[];
    for i=1:1:r(x)
        for a=minbind_hy{x}(i):1:maxbind_hy{x}(i)
            flags=Hydro_1{x}(:,i)>a;
            for f=1:50
                Fl=flags*1;
                Fl(clim(:,f),:)=[];
                Fl_1(:,f)=Fl;
            end
            averghc=(sum(flags(1:(hc),:),1)/(hc)*100);
            averglc=(sum(flags((hc+1):end,:),1)/(lc)*100);
            stdhc=std(sum(Fl_1(1:(hc-nhc)),1)/(hc-nhc)*100);
            stdlc=std(sum(Fl_1((hc-nhc+1):end),1)/(lc-nlc)*100);
            if averghc-stdhc<high && averglc+stdlc>low
               mAb=sum(flags)/(hc+lc)*100;
               highspec=averghc;
               lowspec=averglc;
               flag=a;
               Flags=flags';
               table(1,1)=sum(flags(1:hc)<1);
               table(2,1)=hc-table(1,1);
               table(1,2)=sum(flags(hc+1:end)<1);
               table(2,2)=lc-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=hc+lc;
               pvalue=TwoTailed(table);
               
               Adflags=AdHydro_1{x}(:,i)>a;
               Pflags=PHydro_1{x}(:,i)>a;
               Aflags=AHydro_1{x}(:,i)>a;
               Cflags=CHydro_1{x}(:,i)>a;
               Eflags=EHydro_1{x}(:,i)>a;
               Bflags=BHydro_1{x}(:,i)>a;
               
               table(1,1)=sum(Adflags(1:ha)<1);
               table(2,1)=ha-table(1,1);
               table(1,2)=sum(Adflags(ha+1:end)<1);
               table(2,2)=la-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=ha+la;
               pvalue1=TwoTailed(table);
               
               
               for f=1:50
                   Ad=Adflags;
                   Pflags1=Pflags;
                   Aflags1=Aflags;
                   Cflags1=Cflags;
                   Eflags1=Eflags;
                   Bflags1=Bflags;
                   
                   Ad(adi(:,f),:)=[];
                   Pflags1(clis_p(:,f),:)=[];
                   Aflags1(clis_a(:,f),:)=[];
                   Cflags1(clis_c(:,f),:)=[];
                   Eflags1(clis_e(:,f),:)=[];
                   Bflags1(clis_b(:,f),:)=[];
                   
                   Ad_1(:,f)=Ad;
                   Pflags1_1(:,f)=Pflags1;
                   Aflags1_1(:,f)=Aflags1;
                   Cflags1_1(:,f)=Cflags1;
                   Eflags1_1(:,f)=Eflags1;
                   Bflags1_1(:,f)=Bflags1;
               end
               
               AdflagH=sum(Ad_1(1:(ha-nha),:),1)/(ha-nha)*100;
               AdflagL=sum(Ad_1((ha-nha+1):end,:),1)/(la-nla)*100;
               
               Pflag1H=sum(Pflags1_1(1:(gp-ngp),:),1)/(gp-ngp)*100;
               Aflag1H=sum(Aflags1_1(1:(ga-nga),:),1)/(ga-nga)*100;
               Cflag1H=sum(Cflags1_1(1:(gc-ngc),:),1)/(gc-ngc)*100;
               Eflag1H=sum(Eflags1_1(1:(ge-nge),:),1)/(ge-nge)*100;
               Bflag1H=sum(Bflags1_1(1:(gb-ngb),:),1)/(gb-ngb)*100;
               
               Pflag1L=sum(Pflags1_1((gp-ngp+1):end,:),1)/(bp-nbp)*100;
               Aflag1L=sum(Aflags1_1((ga-nga+1):end,:),1)/(ba-nba)*100;
               Cflag1L=sum(Cflags1_1((gc-ngc+1):end,:),1)/(bc-nbc)*100;
               Eflag1L=sum(Eflags1_1((ge-nge+1):end,:),1)/(be-nbe)*100;
               Bflag1L=sum(Bflags1_1((gb-ngb+1):end,:),1)/(bb-nbb)*100;
               
               avergha=sum(Adflags(1:ha,:),1)/ha*100;
               avergla=sum(Adflags(ha+1:end,:),1)/la*100;
               averghcp=sum(Pflags(1:gp,:),1)/(gp)*100;
               averghca=sum(Aflags(1:ga,:),1)/(ga)*100;
               averghcc=sum(Cflags(1:gc,:),1)/(gc)*100;
               averghce=sum(Eflags(1:ge,:),1)/(ge)*100;
               averghcb=sum(Bflags(1:gb,:),1)/(gb)*100;
               averglcp=sum(Pflags(gp+1:end,:),1)/(bp)*100;
               averglca=sum(Aflags(ga+1:end,:),1)/(ba)*100;
               averglcc=sum(Cflags(gc+1:end,:),1)/(bc)*100;
               averglce=sum(Eflags(ge+1:end,:),1)/(be)*100;
               averglcb=sum(Bflags(gb+1:end,:),1)/(bb)*100;
              
               stdha=std(AdflagH);
               stdla=std(AdflagL);
               stdhcp=std(Pflag1H);
               stdlcp=std(Pflag1L);
               stdhca=std(Aflag1H);
               stdlca=std(Aflag1L);
               stdhcc=std(Cflag1H);
               stdlcc=std(Cflag1L);
               stdhce=std(Eflag1H);
               stdlce=std(Eflag1L);
               stdhcb=std(Bflag1H);
               stdlcb=std(Bflag1L);

                if pvalue<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10 
                  opt{x}(q1,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
                  opt_0=zeros(1,13-x);
                  opt_1{x}(q1,:)=[510,opt_0,Hydro_num(Hydro1{x}(i,:))];
                  opt_2{x}(q1,:)=Flags;
                  q1=q1+1;
                  
               elseif pvalue1<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10  
                  opt{x}(q1,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
                  opt_0=zeros(1,13-x);
                  opt_1{x}(q1,:)=[510,opt_0,Hydro_num(Hydro1{x}(i,:))];
                  opt_2{x}(q1,:)=Flags;
                  q1=q1+1;
               end
            end
        end
    end
end
%% combination 
% polar residue GAVLIMFWP
Polar=[G,A,S,T,Q,N,Y,D,E];
PPolar=[G1,A1,S1,T1,Q1,N1,Y1,D1,E1];
APolar=[G2,A2,S2,T2,Q2,N2,Y2,D2,E2];
CPolar=[G3,A3,S3,T3,Q3,N3,Y3,D3,E3];
EPolar=[G4,A4,S4,T4,Q4,N4,Y4,D4,E4];
BPolar=[G5,A5,S5,T5,Q5,N5,Y5,D5,E5];
AdPolar=[ABG,ABA,ABS,ABT,ABQ,ABN,ABY,ABD,ABE];
Polar_num=[1,2,3,4,5,6,7,8,9];
for i=1:9
    Polar1{i}=combntns(1:9,i);
    [r_p(i),~]=size(Polar1{i});
end 
for k=1:9
    for i=1:1:r_p(k)
        Polar_1{k}(:,i)=sum(Polar(:,Polar1{k}(i,:)),2);
        PPolar_1{k}(:,i)=sum(PPolar(:,Polar1{k}(i,:)),2);
        APolar_1{k}(:,i)=sum(APolar(:,Polar1{k}(i,:)),2);
        CPolar_1{k}(:,i)=sum(CPolar(:,Polar1{k}(i,:)),2);
        EPolar_1{k}(:,i)=sum(EPolar(:,Polar1{k}(i,:)),2);
        BPolar_1{k}(:,i)=sum(BPolar(:,Polar1{k}(i,:)),2);
        AdPolar_1{k}(:,i)=sum(AdPolar(:,Polar1{k}(i,:)),2);
    end
    [~,r_p1(k)]=size(Polar_1{k});
    for j=1:1:r_p1(k)
        maxbind_po{k}(j)=max(Polar_1{k}(:,j));
        minbind_po{k}(j)=min(Polar_1{k}(:,j));
        po{k}(j)=(maxbind_po{k}(j)-minbind_po{k}(j))+1;
    end        
end

for x=1:1:9
    q2=1;
    opt_0=[];
    for i=1:1:r_p(x)
        for a=minbind_po{x}(i):1:maxbind_po{x}(i)
             flags=Polar_1{x}(:,i)<a;
            for f=1:50
                Fl=flags*1;
                Fl(clim(:,f),:)=[];
                Fl_1(:,f)=Fl;
            end
            averghc=sum(flags(1:(hc),:),1)/(hc)*100;
            averglc=(sum(flags((hc+1):end,:),1)/(lc)*100);
            stdhc=std(sum(Fl_1(1:(hc-nhc)),1)/(hc-nhc)*100);
            stdlc=std(sum(Fl_1((hc-nhc+1):end),1)/(lc-nlc)*100);
            if averghc-stdhc<high && averglc+stdlc>low
               mAb=sum(flags)/(hc+lc)*100;
               highspec=averghc;
               lowspec=averglc;
               flag=a;
               Flags=flags';
               table(1,1)=sum(flags(1:hc)<1);
               table(2,1)=hc-table(1,1);
               table(1,2)=sum(flags(hc+1:end)<1);
               table(2,2)=lc-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=hc+lc;
               pvalue=TwoTailed(table);
               
               Adflags=AdPolar_1{x}(:,i)<a;
               Pflags=PPolar_1{x}(:,i)<a;
               Aflags=APolar_1{x}(:,i)<a;
               Cflags=CPolar_1{x}(:,i)<a;
               Eflags=EPolar_1{x}(:,i)<a;
               Bflags=BPolar_1{x}(:,i)<a;
               
               table(1,1)=sum(Adflags(1:ha)<1);
               table(2,1)=ha-table(1,1);
               table(1,2)=sum(Adflags(ha+1:end)<1);
               table(2,2)=la-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=ha+la;
               pvalue1=TwoTailed(table);
               
               for f=1:50
                   Ad=Adflags;
                   Pflags1=Pflags;
                   Aflags1=Aflags;
                   Cflags1=Cflags;
                   Eflags1=Eflags;
                   Bflags1=Bflags;
                   
                   Ad(adi(:,f),:)=[];
                   Pflags1(clis_p(:,f),:)=[];
                   Aflags1(clis_a(:,f),:)=[];
                   Cflags1(clis_c(:,f),:)=[];
                   Eflags1(clis_e(:,f),:)=[];
                   Bflags1(clis_b(:,f),:)=[];
                   
                   Ad_1(:,f)=Ad;
                   Pflags1_1(:,f)=Pflags1;
                   Aflags1_1(:,f)=Aflags1;
                   Cflags1_1(:,f)=Cflags1;
                   Eflags1_1(:,f)=Eflags1;
                   Bflags1_1(:,f)=Bflags1;
               end
               
               AdflagH=sum(Ad_1(1:(ha-nha),:),1)/(ha-nha)*100;
               AdflagL=sum(Ad_1((ha-nha+1):end,:),1)/(la-nla)*100;
               
               Pflag1H=sum(Pflags1_1(1:(gp-ngp),:),1)/(gp-ngp)*100;
               Aflag1H=sum(Aflags1_1(1:(ga-nga),:),1)/(ga-nga)*100;
               Cflag1H=sum(Cflags1_1(1:(gc-ngc),:),1)/(gc-ngc)*100;
               Eflag1H=sum(Eflags1_1(1:(ge-nge),:),1)/(ge-nge)*100;
               Bflag1H=sum(Bflags1_1(1:(gb-ngb),:),1)/(gb-ngb)*100;
               
               Pflag1L=sum(Pflags1_1((gp-ngp+1):end,:),1)/(bp-nbp)*100;
               Aflag1L=sum(Aflags1_1((ga-nga+1):end,:),1)/(ba-nba)*100;
               Cflag1L=sum(Cflags1_1((gc-ngc+1):end,:),1)/(bc-nbc)*100;
               Eflag1L=sum(Eflags1_1((ge-nge+1):end,:),1)/(be-nbe)*100;
               Bflag1L=sum(Bflags1_1((gb-ngb+1):end,:),1)/(bb-nbb)*100;
               
               avergha=sum(Adflags(1:ha,:),1)/ha*100;
               avergla=sum(Adflags(ha+1:end,:),1)/la*100;
               averghcp=sum(Pflags(1:gp,:),1)/(gp)*100;
               averghca=sum(Aflags(1:ga,:),1)/(ga)*100;
               averghcc=sum(Cflags(1:gc,:),1)/(gc)*100;
               averghce=sum(Eflags(1:ge,:),1)/(ge)*100;
               averghcb=sum(Bflags(1:gb,:),1)/(gb)*100;
               averglcp=sum(Pflags(gp+1:end,:),1)/(bp)*100;
               averglca=sum(Aflags(ga+1:end,:),1)/(ba)*100;
               averglcc=sum(Cflags(gc+1:end,:),1)/(bc)*100;
               averglce=sum(Eflags(ge+1:end,:),1)/(be)*100;
               averglcb=sum(Bflags(gb+1:end,:),1)/(bb)*100;
              
               stdha=std(AdflagH);
               stdla=std(AdflagL);
               stdhcp=std(Pflag1H);
               stdlcp=std(Pflag1L);
               stdhca=std(Aflag1H);
               stdlca=std(Aflag1L);
               stdhcc=std(Cflag1H);
               stdlcc=std(Cflag1L);
               stdhce=std(Eflag1H);
               stdlce=std(Eflag1L);
               stdhcb=std(Bflag1H);
               stdlcb=std(Bflag1L);  
               if pvalue<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10 
                  opt{x+13}(q2,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
                  opt_0=zeros(1,13-x);
                  opt_1{x+13}(q2,:)=[100,opt_0,Polar_num(Polar1{x}(i,:))];
                  opt_2{x+13}(q2,:)=Flags;
                  q2=q2+1;
               elseif pvalue1<0.05 &&averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10  
                  opt{x+13}(q2,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
                  opt_0=zeros(1,13-x);
                  opt_1{x+13}(q2,:)=[100,opt_0,Polar_num(Polar1{x}(i,:))];
                  opt_2{x+13}(q2,:)=Flags;
                  q2=q2+1;
               end
            end
        end
    end
end
%% single amino acid
Single=netcharge;
PSingle=Pnetcharge;
ASingle=Anetcharge;
CSingle=Cnetcharge;
ESingle=Enetcharge;
BSingle=Bnetcharge;
AdSingle=ABnetcharge;
y=1;
for a=0:1:max_nc
    flags=Single>a;
    for f=1:50
        Fl=flags*1;
        Fl(clim(:,f),:)=[];
        Fl_1(:,f)=Fl;
    end
    averghc=(sum(flags(1:(hc),:),1)/(hc)*100);
            averglc=(sum(flags((hc+1):end,:),1)/(lc)*100);
            stdhc=std(sum(Fl_1(1:(hc-nhc)),1)/(hc-nhc)*100);
            stdlc=std(sum(Fl_1((hc-nhc+1):end),1)/(lc-nlc)*100);
            if averghc-stdhc<high && averglc+stdlc>low
               mAb=sum(flags)/(hc+lc)*100;
               highspec=averghc;
               lowspec=averglc;
               flag=a;
               Flags=flags';
               table(1,1)=sum(flags(1:hc)<1);
               table(2,1)=hc-table(1,1);
               table(1,2)=sum(flags(hc+1:end)<1);
               table(2,2)=lc-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=hc+lc;
               pvalue=TwoTailed(table);
       
       Adflags=AdSingle>a;
       Pflags=PSingle>a;
       Aflags=ASingle>a;
       Cflags=CSingle>a;
       Eflags=ESingle>a;
       Bflags=BSingle>a;
       
       table(1,1)=sum(Adflags(1:ha)<1);
       table(2,1)=ha-table(1,1);
       table(1,2)=sum(Adflags(ha+1:end)<1);
       table(2,2)=la-table(1,2);
       table(3,1)=table(1,1)+table(2,1);
       table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=ha+la;
               pvalue1=TwoTailed(table);
       for f=1:50
                   Ad=Adflags;
                   Pflags1=Pflags;
                   Aflags1=Aflags;
                   Cflags1=Cflags;
                   Eflags1=Eflags;
                   Bflags1=Bflags;
                   
                   Ad(adi(:,f),:)=[];
                   Pflags1(clis_p(:,f),:)=[];
                   Aflags1(clis_a(:,f),:)=[];
                   Cflags1(clis_c(:,f),:)=[];
                   Eflags1(clis_e(:,f),:)=[];
                   Bflags1(clis_b(:,f),:)=[];
                   
                   Ad_1(:,f)=Ad;
                   Pflags1_1(:,f)=Pflags1;
                   Aflags1_1(:,f)=Aflags1;
                   Cflags1_1(:,f)=Cflags1;
                   Eflags1_1(:,f)=Eflags1;
                   Bflags1_1(:,f)=Bflags1;
               end
               
               AdflagH=sum(Ad_1(1:(ha-nha),:),1)/(ha-nha)*100;
               AdflagL=sum(Ad_1((ha-nha+1):end,:),1)/(la-nla)*100;
               
               Pflag1H=sum(Pflags1_1(1:(gp-ngp),:),1)/(gp-ngp)*100;
               Aflag1H=sum(Aflags1_1(1:(ga-nga),:),1)/(ga-nga)*100;
               Cflag1H=sum(Cflags1_1(1:(gc-ngc),:),1)/(gc-ngc)*100;
               Eflag1H=sum(Eflags1_1(1:(ge-nge),:),1)/(ge-nge)*100;
               Bflag1H=sum(Bflags1_1(1:(gb-ngb),:),1)/(gb-ngb)*100;
               
               Pflag1L=sum(Pflags1_1((gp-ngp+1):end,:),1)/(bp-nbp)*100;
               Aflag1L=sum(Aflags1_1((ga-nga+1):end,:),1)/(ba-nba)*100;
               Cflag1L=sum(Cflags1_1((gc-ngc+1):end,:),1)/(bc-nbc)*100;
               Eflag1L=sum(Eflags1_1((ge-nge+1):end,:),1)/(be-nbe)*100;
               Bflag1L=sum(Bflags1_1((gb-ngb+1):end,:),1)/(bb-nbb)*100;
               
               avergha=sum(Adflags(1:ha,:),1)/ha*100;
               avergla=sum(Adflags(ha+1:end,:),1)/la*100;
               averghcp=sum(Pflags(1:gp,:),1)/(gp)*100;
               averghca=sum(Aflags(1:ga,:),1)/(ga)*100;
               averghcc=sum(Cflags(1:gc,:),1)/(gc)*100;
               averghce=sum(Eflags(1:ge,:),1)/(ge)*100;
               averghcb=sum(Bflags(1:gb,:),1)/(gb)*100;
               averglcp=sum(Pflags(gp+1:end,:),1)/(bp)*100;
               averglca=sum(Aflags(ga+1:end,:),1)/(ba)*100;
               averglcc=sum(Cflags(gc+1:end,:),1)/(bc)*100;
               averglce=sum(Eflags(ge+1:end,:),1)/(be)*100;
               averglcb=sum(Bflags(gb+1:end,:),1)/(bb)*100;
              
               stdha=std(AdflagH);
               stdla=std(AdflagL);
               stdhcp=std(Pflag1H);
               stdlcp=std(Pflag1L);
               stdhca=std(Aflag1H);
               stdlca=std(Aflag1L);
               stdhcc=std(Cflag1H);
               stdlcc=std(Cflag1L);
               stdhce=std(Eflag1H);
               stdlce=std(Eflag1L);
               stdhcb=std(Bflag1H);
               stdlcb=std(Bflag1L);
       if pvalue<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10  
          opt{45}(y,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
          opt_1{45}(y,:)=[710,0,0,0,0,0,0,0,0,0,0,0,0,1];
          opt_2{45}(y,:)=[Flags];
          y=y+1;
       elseif pvalue1<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10  
          opt{45}(y,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
          opt_1{45}(y,:)=[710,0,0,0,0,0,0,0,0,0,0,0,0,1];
          opt_2{45}(y,:)=[Flags];
          y=y+1;
       end
    end
end
Row=[];               
exist opt               
if ans~=0
   [~,r]=size(opt);
   Res=[];
   for i=1:r
       Res=[Res;opt_1{i} opt{i} opt_2{i}];
   end
   Res(isnan(Res(:,1)),:) = [] ;

   [e,~]=size(Res);
   Fx=Res(:,1:14);
   Px=Res(:,15:19);
   q=0;
   for i=1:e-1
       for j=i+1:e
           p=ismember(Fx(i,:),Fx(j,:));
           p1=sum(p);
           if p1==14 && Px(i,1)==Px(j,1)&& Px(i,2)==Px(j,2) && Px(i,3)==Px(j,3)
              q=q+1;
              Row(q)=j;
           end
       end
   end

   if q~=0
      Row=unique(Row);
      Res(Row,:) = []; 
      Res(:,18)=-1*Res(:,18);
      Res(:,20)=-1*Res(:,20);
      Res(:,22)=-1*Res(:,22);
      Res(:,24)=-1*Res(:,24);
      Res(:,26)=-1*Res(:,26);
      Res(:,28)=-1*Res(:,28);
      Res(:,30)=-1*Res(:,30);
      Res=sortrows(Res,[16,18,30,20,22,24,26,28]);
      Res(:,18)=-1*Res(:,18);
      Res(:,20)=-1*Res(:,20);
      Res(:,22)=-1*Res(:,22);
      Res(:,24)=-1*Res(:,24);
      Res(:,26)=-1*Res(:,26);
      Res(:,28)=-1*Res(:,28);
      Res(:,30)=-1*Res(:,30);
      [res,~]=size(Res);
      Sig=zeros(res,1);
      for i=1:res
          pc=i/res*0.05;
          if Res(i,16)<pc
             Sig(i,:)=1;
          end
      end
      s=find(Sig==1);
      [rs,~]=size(s);
      if rs==0
         disp('No results')
      else
      max1=max(s);
      Res_1=Res(1:max1,:); % results from FDR
      end
   else
      Res(:,18)=-1*Res(:,18);
      Res(:,20)=-1*Res(:,20);
      Res(:,22)=-1*Res(:,22);
      Res(:,24)=-1*Res(:,24);
      Res(:,26)=-1*Res(:,26);
      Res(:,28)=-1*Res(:,28);
      Res(:,30)=-1*Res(:,30);
      Res=sortrows(Res,[16,18,30,20,22,24,26,28]);
      Res(:,18)=-1*Res(:,18);
      Res(:,20)=-1*Res(:,20);
      Res(:,22)=-1*Res(:,22);
      Res(:,24)=-1*Res(:,24);
      Res(:,26)=-1*Res(:,26);
      Res(:,28)=-1*Res(:,28);
      Res(:,30)=-1*Res(:,30);
      [res,~]=size(Res);
      Sig=zeros(res,1);
      for i=1:res
          pc=i/res*0.05;
          if Res(i,16)<pc
          Sig(i,:)=1;
          end
      end
      s=find(Sig==1);
      [rs,~]=size(s);
      if rs==0
         disp('No results')
      else
         max1=max(s);
         Res_1=Res(1:max1,:); % results from FDR
      end
   end
   exist Res_1
   if ans ~=0 
      [res1,~]=size(Res_1);
      if res1~=0
         Sign_v=Res_1(:,1);
         Flag_AA=Res_1(:,2:14);
         other_info=Res_1(:,15:end);
         for i=1:res1
             F=nonzeros(Flag_AA(i,:))';
             if Sign_v(i)==510
                Hydro_AA='GAVLIMFWPYRKH';
                Flag_AA1{i}=Hydro_AA(F);
                Sign1{i}='>';
             elseif Sign_v(i)==100
                Pola_AA='GASTQNYDE';
                Flag_AA1{i}=Pola_AA(F);
                Sign1{i}='<';
             elseif Sign_v(i)==710
                Flag_AA1{i}='charge';
                Sign1{i}='>';
             end
         end
      end
      Flag=cell2table(Flag_AA1');
      Sign=cell2table(Sign1');
      Flag.Properties.VariableNames{'Var1'}='Flag';
      Sign.Properties.VariableNames{'Var1'}='Sign';
      info=array2table(other_info);
      R=cell(res1,1);
      R(:)={Region};
      Singleflag=[R,Flag,Sign,info];
   else
      Singleflag=[];
   end
   
  else
   Singleflag=[];
end
region
  Finalflag=[Finalflag;Singleflag];
end
Region1={'H1','H1','H1','H1','H1','H2','H2','H2','H2','H3','H3','H3','L1','L1','L2'};
Region2={'H2','H3','L1','L2','L3','H3','L1','L2','L3','L1','L2','L3','L2','L3','L3'};
Region={'H1H2','H1H3','H1L1','H1L2','H1L3','H2H3','H2L1','H2L2','H2L3','H3L1','H3L2','H3L3','L1L2','L1L3','L2L3'};
for re2=1:15
    clear opt opt_0 opt_1 opt_2 Flag_AA1 Sign1 Sign_v Res_1 Row
warning('off','all');
CAB1=xlsread('Clinical_AA',Region1{re2},'A1:T200');
%HAB=xlsread('human antibody_sequence_AA comp_newlimit_v1',Region,'A2:T20013');
CP1=xlsread('Clinical_AA_PSR',Region1{re2},'A1:T300');
CA1=xlsread('Clinical_AA_ACSINS',Region1{re2},'A1:T300');
CC1=xlsread('Clinical_AA_CSI',Region1{re2},'A1:T300');
CE1=xlsread('Clinical_AA_ELISA',Region1{re2},'A1:T300');
CB1=xlsread('Clinical_AA_BVP',Region1{re2},'A1:T300');
AAB1=xlsread('Human_AA_100sim',Region1{re2},'A1:T2000');

CAB2=xlsread('Clinical_AA',Region2{re2},'A1:T200');
%HAB=xlsread('human antibody_sequence_AA comp_newlimit_v1',Region,'A2:T20013');
CP2=xlsread('Clinical_AA_PSR',Region2{re2},'A1:T300');
CA2=xlsread('Clinical_AA_ACSINS',Region2{re2},'A1:T300');
CC2=xlsread('Clinical_AA_CSI',Region2{re2},'A1:T300');
CE2=xlsread('Clinical_AA_ELISA',Region2{re2},'A1:T300');
CB2=xlsread('Clinical_AA_BVP',Region2{re2},'A1:T300');
AAB2=xlsread('Human_AA_100sim',Region2{re2},'A1:T2000');

CAB=CAB1+CAB2;
CP=CP1+CP2;
CA=CA1+CA2;
CC=CC1+CC2;
CE=CE1+CE2;
CB=CB1+CB2;
AAB=AAB1+AAB2;
%hv=1033;lv=59;hr=104;lr=6;highP=0.15;
adihs=[];adils=[];clis_ph=[];clis_pl=[];clis_ch=[];clis_cl=[];clis_eh=[];clis_el=[];clis_bh=[];clis_bl=[];
clis_ah=[];clis_al=[];clihs=[];clils=[];

adi=xlsread('2ndset5_combine_singleflags_v2_new','adimab','A2:CV100');
clis_p=xlsread('2ndset5_combine_singleflags_v2_new','PSR','A2:CV20');
clis_a=xlsread('2ndset5_combine_singleflags_v2_new','ACSINS','A2:CV20');
clis_c=xlsread('2ndset5_combine_singleflags_v2_new','CSI','A2:CV20');
clis_e=xlsread('2ndset5_combine_singleflags_v2_new','ELISA','A2:CV20');
clis_b=xlsread('2ndset5_combine_singleflags_v2_new','BVP','A2:CV20');
clim=xlsread('2ndset5_combine_singleflags_v2_new','clinical','A2:CV100');
ngp=sum(clis_p(:,1)<=gp);nbp=sum(clis_p(:,1)>gp);
nga=sum(clis_a(:,1)<=ga);nba=sum(clis_a(:,1)>ga);
ngc=sum(clis_c(:,1)<=gc);nbc=sum(clis_c(:,1)>gc);
nge=sum(clis_e(:,1)<=ge);nbe=sum(clis_e(:,1)>ge);
ngb=sum(clis_b(:,1)<=gb);nbb=sum(clis_b(:,1)>gb);
nha=sum(adi(:,1)<=ha);nla=sum(adi(:,1)>ha);
nhc=sum(clim(:,1)<=hc);nlc=sum(clim(:,1)>hc);
%hv=1033;lv=59;hr=104;lr=6;highP=0.15;
A=CAB(Total2,1);C=CAB(Total2,2);D=CAB(Total2,3);E=CAB(Total2,4);F=CAB(Total2,5);G=CAB(Total2,6);H=CAB(Total2,7);
I=CAB(Total2,8);K=CAB(Total2,9);L=CAB(Total2,10);M=CAB(Total2,11);N=CAB(Total2,12);P=CAB(Total2,13);Q=CAB(Total2,14);
R=CAB(Total2,15);S=CAB(Total2,16);T=CAB(Total2,17);V=CAB(Total2,18);W=CAB(Total2,19);Y=CAB(Total2,20);
netcharge=R+K+0.1*H-D-E;
max_nc=floor(max(netcharge));
%% Human Antibody information PSR
A1=CP(PSR2,1);C1=CP(PSR2,2);D1=CP(PSR2,3);E1=CP(PSR2,4);F1=CP(PSR2,5);G1=CP(PSR2,6);H1=CP(PSR2,7);
I1=CP(PSR2,8);K1=CP(PSR2,9);L1=CP(PSR2,10);M1=CP(PSR2,11);N1=CP(PSR2,12);P1=CP(PSR2,13);Q1=CP(PSR2,14);
R1=CP(PSR2,15);S1=CP(PSR2,16);T1=CP(PSR2,17);V1=CP(PSR2,18);W1=CP(PSR2,19);Y1=CP(PSR2,20);
Pnetcharge=R1+K1+0.1*H1-D1-E1;
%% Human Antibody information ACSINS
A2=CA(ACSINS2,1);C2=CA(ACSINS2,2);D2=CA(ACSINS2,3);E2=CA(ACSINS2,4);F2=CA(ACSINS2,5);G2=CA(ACSINS2,6);H2=CA(ACSINS2,7);
I2=CA(ACSINS2,8);K2=CA(ACSINS2,9);L2=CA(ACSINS2,10);M2=CA(ACSINS2,11);N2=CA(ACSINS2,12);P2=CA(ACSINS2,13);Q2=CA(ACSINS2,14);
R2=CA(ACSINS2,15);S2=CA(ACSINS2,16);T2=CA(ACSINS2,17);V2=CA(ACSINS2,18);W2=CA(ACSINS2,19);Y2=CA(ACSINS2,20);
Anetcharge=R2+K2+0.1*H2-D2-E2;
%% Human Antibody information CSI
A3=CC(CSI2,1);D3=CC(CSI2,3);E3=CC(CSI2,4);F3=CC(CSI2,5);G3=CC(CSI2,6);H3=CC(CSI2,7);
I3=CC(CSI2,8);K3=CC(CSI2,9);L3=CC(CSI2,10);M3=CC(CSI2,11);N3=CC(CSI2,12);P3=CC(CSI2,13);Q3=CC(CSI2,14);
R3=CC(CSI2,15);S3=CC(CSI2,16);T3=CC(CSI2,17);V3=CC(CSI2,18);W3=CC(CSI2,19);Y3=CC(CSI2,20);
Cnetcharge=R3+K3+0.1*H3-D3-E3;
%% Human Antibody information ELISA
A4=CE(ELISA2,1);C4=CE(ELISA2,2);D4=CE(ELISA2,3);E4=CE(ELISA2,4);F4=CE(ELISA2,5);G4=CE(ELISA2,6);H4=CE(ELISA2,7);
I4=CE(ELISA2,8);K4=CE(ELISA2,9);L4=CE(ELISA2,10);M4=CE(ELISA2,11);N4=CE(ELISA2,12);P4=CE(ELISA2,13);Q4=CE(ELISA2,14);
R4=CE(ELISA2,15);S4=CE(ELISA2,16);T4=CE(ELISA2,17);V4=CE(ELISA2,18);W4=CE(ELISA2,19);Y4=CE(ELISA2,20);
Enetcharge=R4+K4+0.1*H4-D4-E4;
%% Human Antibody information ELISA
A5=CB(BVP2,1);C5=CB(BVP2,2);D5=CB(BVP2,3);E5=CB(BVP2,4);F5=CB(BVP2,5);G5=CB(BVP2,6);H5=CB(BVP2,7);
I5=CB(BVP2,8);K5=CB(BVP2,9);L5=CB(BVP2,10);M5=CB(BVP2,11);N5=CB(BVP2,12);P5=CB(BVP2,13);Q5=CB(BVP2,14);
R5=CB(BVP2,15);S5=CB(BVP2,16);T5=CB(BVP2,17);V5=CB(BVP2,18);W5=CB(BVP2,19);Y5=CB(BVP2,20);
Bnetcharge=R5+K5+0.1*H5-D5-E5;
%% Human B cell from Adimab
ABA=AAB(train2,1);ABC=AAB(train2,2);ABD=AAB(train2,3);ABE=AAB(train2,4);ABF=AAB(train2,5);ABG=AAB(train2,6);ABH=AAB(train2,7);
ABI=AAB(train2,8);ABK=AAB(train2,9);ABL=AAB(train2,10);ABM=AAB(train2,11);ABN=AAB(train2,12);ABP=AAB(train2,13);ABQ=AAB(train2,14);
ABR=AAB(train2,15);ABS=AAB(train2,16);ABT=AAB(train2,17);ABV=AAB(train2,18);ABW=AAB(train2,19);ABY=AAB(train2,20);
ABnetcharge=ABR+ABK+0.1*ABH-ABD-ABE;
%% combination
% hydrophybic residue GAVLIMFWPYRKH
Hydro=[G,A,V,L,I,M,F,W,P,Y,R,K,H];
PHydro=[G1,A1,V1,L1,I1,M1,F1,W1,P1,Y1,R1,K1,H1];
AHydro=[G2,A2,V2,L2,I2,M2,F2,W2,P2,Y2,R2,K2,H2];
CHydro=[G3,A3,V3,L3,I3,M3,F3,W3,P3,Y3,R3,K3,H3];
EHydro=[G4,A4,V4,L4,I4,M4,F4,W4,P4,Y4,R4,K4,H4];
BHydro=[G5,A5,V5,L5,I5,M5,F5,W5,P5,Y5,R5,K5,H5];
AdHydro=[ABG,ABA,ABV,ABL,ABI,ABM,ABF,ABW,ABP,ABY,ABR,ABK,ABH];
Hydro_num=[1,2,3,4,5,6,7,8,9,10,11,12,13];
for i=1:13
    Hydro1{i}=combntns(1:13,i);
    [r(i),~]=size(Hydro1{i});
end    
for k=1:13
    for i=1:1:r(k)
        Hydro_1{k}(:,i)=sum(Hydro(:,Hydro1{k}(i,:)),2);
        PHydro_1{k}(:,i)=sum(PHydro(:,Hydro1{k}(i,:)),2);
        AHydro_1{k}(:,i)=sum(AHydro(:,Hydro1{k}(i,:)),2);
        CHydro_1{k}(:,i)=sum(CHydro(:,Hydro1{k}(i,:)),2);
        EHydro_1{k}(:,i)=sum(EHydro(:,Hydro1{k}(i,:)),2);
        BHydro_1{k}(:,i)=sum(BHydro(:,Hydro1{k}(i,:)),2);
        AdHydro_1{k}(:,i)=sum(AdHydro(:,Hydro1{k}(i,:)),2);
    end
    [~,r_1(k)]=size(Hydro_1{k});
    for j=1:1:r_1(k)
        maxbind_hy{k}(j)=max(Hydro_1{k}(:,j));
        minbind_hy{k}(j)=min(Hydro_1{k}(:,j));
        hy{k}(j)=(maxbind_hy{k}(j)-minbind_hy{k}(j))+1;
    end        
end    
for x=1:1:13
    q1=1;
    opt_0=[];
    for i=1:1:r(x)
        for a=minbind_hy{x}(i):1:maxbind_hy{x}(i)
            flags=Hydro_1{x}(:,i)>a;
            for f=1:50
                Fl=flags*1;
                Fl(clim(:,f),:)=[];
                Fl_1(:,f)=Fl;
            end
            averghc=(sum(flags(1:(hc),:),1)/(hc)*100);
            averglc=(sum(flags((hc+1):end,:),1)/(lc)*100);
            stdhc=std(sum(Fl_1(1:(hc-nhc)),1)/(hc-nhc)*100);
            stdlc=std(sum(Fl_1((hc-nhc+1):end),1)/(lc-nlc)*100);
            if averghc-stdhc<high && averglc+stdlc>low
               mAb=sum(flags)/(hc+lc)*100;
               highspec=averghc;
               lowspec=averglc;
               flag=a;
               Flags=flags';
               table(1,1)=sum(flags(1:hc)<1);
               table(2,1)=hc-table(1,1);
               table(1,2)=sum(flags(hc+1:end)<1);
               table(2,2)=lc-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=hc+lc;
               pvalue=TwoTailed(table);
               
               Adflags=AdHydro_1{x}(:,i)>a;
               Pflags=PHydro_1{x}(:,i)>a;
               Aflags=AHydro_1{x}(:,i)>a;
               Cflags=CHydro_1{x}(:,i)>a;
               Eflags=EHydro_1{x}(:,i)>a;
               Bflags=BHydro_1{x}(:,i)>a;
               
               table(1,1)=sum(Adflags(1:ha)<1);
               table(2,1)=ha-table(1,1);
               table(1,2)=sum(Adflags(ha+1:end)<1);
               table(2,2)=la-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=ha+la;
               pvalue1=TwoTailed(table);
               
               
               for f=1:50
                   Ad=Adflags;
                   Pflags1=Pflags;
                   Aflags1=Aflags;
                   Cflags1=Cflags;
                   Eflags1=Eflags;
                   Bflags1=Bflags;
                   
                   Ad(adi(:,f),:)=[];
                   Pflags1(clis_p(:,f),:)=[];
                   Aflags1(clis_a(:,f),:)=[];
                   Cflags1(clis_c(:,f),:)=[];
                   Eflags1(clis_e(:,f),:)=[];
                   Bflags1(clis_b(:,f),:)=[];
                   
                   Ad_1(:,f)=Ad;
                   Pflags1_1(:,f)=Pflags1;
                   Aflags1_1(:,f)=Aflags1;
                   Cflags1_1(:,f)=Cflags1;
                   Eflags1_1(:,f)=Eflags1;
                   Bflags1_1(:,f)=Bflags1;
               end
               
               AdflagH=sum(Ad_1(1:(ha-nha),:),1)/(ha-nha)*100;
               AdflagL=sum(Ad_1((ha-nha+1):end,:),1)/(la-nla)*100;
               
               Pflag1H=sum(Pflags1_1(1:(gp-ngp),:),1)/(gp-ngp)*100;
               Aflag1H=sum(Aflags1_1(1:(ga-nga),:),1)/(ga-nga)*100;
               Cflag1H=sum(Cflags1_1(1:(gc-ngc),:),1)/(gc-ngc)*100;
               Eflag1H=sum(Eflags1_1(1:(ge-nge),:),1)/(ge-nge)*100;
               Bflag1H=sum(Bflags1_1(1:(gb-ngb),:),1)/(gb-ngb)*100;
               
               Pflag1L=sum(Pflags1_1((gp-ngp+1):end,:),1)/(bp-nbp)*100;
               Aflag1L=sum(Aflags1_1((ga-nga+1):end,:),1)/(ba-nba)*100;
               Cflag1L=sum(Cflags1_1((gc-ngc+1):end,:),1)/(bc-nbc)*100;
               Eflag1L=sum(Eflags1_1((ge-nge+1):end,:),1)/(be-nbe)*100;
               Bflag1L=sum(Bflags1_1((gb-ngb+1):end,:),1)/(bb-nbb)*100;
               
               avergha=sum(Adflags(1:ha,:),1)/ha*100;
               avergla=sum(Adflags(ha+1:end,:),1)/la*100;
               averghcp=sum(Pflags(1:gp,:),1)/(gp)*100;
               averghca=sum(Aflags(1:ga,:),1)/(ga)*100;
               averghcc=sum(Cflags(1:gc,:),1)/(gc)*100;
               averghce=sum(Eflags(1:ge,:),1)/(ge)*100;
               averghcb=sum(Bflags(1:gb,:),1)/(gb)*100;
               averglcp=sum(Pflags(gp+1:end,:),1)/(bp)*100;
               averglca=sum(Aflags(ga+1:end,:),1)/(ba)*100;
               averglcc=sum(Cflags(gc+1:end,:),1)/(bc)*100;
               averglce=sum(Eflags(ge+1:end,:),1)/(be)*100;
               averglcb=sum(Bflags(gb+1:end,:),1)/(bb)*100;
              
               stdha=std(AdflagH);
               stdla=std(AdflagL);
               stdhcp=std(Pflag1H);
               stdlcp=std(Pflag1L);
               stdhca=std(Aflag1H);
               stdlca=std(Aflag1L);
               stdhcc=std(Cflag1H);
               stdlcc=std(Cflag1L);
               stdhce=std(Eflag1H);
               stdlce=std(Eflag1L);
               stdhcb=std(Bflag1H);
               stdlcb=std(Bflag1L);
               
                  
               if pvalue<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10 
                  opt{x}(q1,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
                  opt_0=zeros(1,13-x);
                  opt_1{x}(q1,:)=[510,opt_0,Hydro_num(Hydro1{x}(i,:))];
                  opt_2{x}(q1,:)=Flags;
                  q1=q1+1;
                  
               elseif pvalue1<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10  
                  opt{x}(q1,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
                  opt_0=zeros(1,13-x);
                  opt_1{x}(q1,:)=[510,opt_0,Hydro_num(Hydro1{x}(i,:))];
                  opt_2{x}(q1,:)=Flags;
                  q1=q1+1;
               end
            end
        end
    end
end
%% combination 
% polar residue GAVLIMFWP
Polar=[G,A,S,T,Q,N,Y,D,E];
PPolar=[G1,A1,S1,T1,Q1,N1,Y1,D1,E1];
APolar=[G2,A2,S2,T2,Q2,N2,Y2,D2,E2];
CPolar=[G3,A3,S3,T3,Q3,N3,Y3,D3,E3];
EPolar=[G4,A4,S4,T4,Q4,N4,Y4,D4,E4];
BPolar=[G5,A5,S5,T5,Q5,N5,Y5,D5,E5];
AdPolar=[ABG,ABA,ABS,ABT,ABQ,ABN,ABY,ABD,ABE];
Polar_num=[1,2,3,4,5,6,7,8,9];
for i=1:9
    Polar1{i}=combntns(1:9,i);
    [r_p(i),~]=size(Polar1{i});
end 
for k=1:9
    for i=1:1:r_p(k)
        Polar_1{k}(:,i)=sum(Polar(:,Polar1{k}(i,:)),2);
        PPolar_1{k}(:,i)=sum(PPolar(:,Polar1{k}(i,:)),2);
        APolar_1{k}(:,i)=sum(APolar(:,Polar1{k}(i,:)),2);
        CPolar_1{k}(:,i)=sum(CPolar(:,Polar1{k}(i,:)),2);
        EPolar_1{k}(:,i)=sum(EPolar(:,Polar1{k}(i,:)),2);
        BPolar_1{k}(:,i)=sum(BPolar(:,Polar1{k}(i,:)),2);
        AdPolar_1{k}(:,i)=sum(AdPolar(:,Polar1{k}(i,:)),2);
    end
    [~,r_p1(k)]=size(Polar_1{k});
    for j=1:1:r_p1(k)
        maxbind_po{k}(j)=max(Polar_1{k}(:,j));
        minbind_po{k}(j)=min(Polar_1{k}(:,j));
        po{k}(j)=(maxbind_po{k}(j)-minbind_po{k}(j))+1;
    end        
end

for x=1:1:9
    q2=1;
    opt_0=[];
    for i=1:1:r_p(x)
        for a=minbind_po{x}(i):1:maxbind_po{x}(i)
             flags=Polar_1{x}(:,i)<a;
            for f=1:50
                Fl=flags*1;
                Fl(clim(:,f),:)=[];
                Fl_1(:,f)=Fl;
            end
            averghc=sum(flags(1:(hc),:),1)/(hc)*100;
            averglc=(sum(flags((hc+1):end,:),1)/(lc)*100);
            stdhc=std(sum(Fl_1(1:(hc-nhc)),1)/(hc-nhc)*100);
            stdlc=std(sum(Fl_1((hc-nhc+1):end),1)/(lc-nlc)*100);
            if averghc-stdhc<high && averglc+stdlc>low
               mAb=sum(flags)/(hc+lc)*100;
               highspec=averghc;
               lowspec=averglc;
               flag=a;
               Flags=flags';
               table(1,1)=sum(flags(1:hc)<1);
               table(2,1)=hc-table(1,1);
               table(1,2)=sum(flags(hc+1:end)<1);
               table(2,2)=lc-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=hc+lc;
               pvalue=TwoTailed(table);
               
               Adflags=AdPolar_1{x}(:,i)<a;
               Pflags=PPolar_1{x}(:,i)<a;
               Aflags=APolar_1{x}(:,i)<a;
               Cflags=CPolar_1{x}(:,i)<a;
               Eflags=EPolar_1{x}(:,i)<a;
               Bflags=BPolar_1{x}(:,i)<a;
               
               table(1,1)=sum(Adflags(1:ha)<1);
               table(2,1)=ha-table(1,1);
               table(1,2)=sum(Adflags(ha+1:end)<1);
               table(2,2)=la-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=ha+la;
               pvalue1=TwoTailed(table);
               
               for f=1:50
                   Ad=Adflags;
                   Pflags1=Pflags;
                   Aflags1=Aflags;
                   Cflags1=Cflags;
                   Eflags1=Eflags;
                   Bflags1=Bflags;
                   
                   Ad(adi(:,f),:)=[];
                   Pflags1(clis_p(:,f),:)=[];
                   Aflags1(clis_a(:,f),:)=[];
                   Cflags1(clis_c(:,f),:)=[];
                   Eflags1(clis_e(:,f),:)=[];
                   Bflags1(clis_b(:,f),:)=[];
                   
                   Ad_1(:,f)=Ad;
                   Pflags1_1(:,f)=Pflags1;
                   Aflags1_1(:,f)=Aflags1;
                   Cflags1_1(:,f)=Cflags1;
                   Eflags1_1(:,f)=Eflags1;
                   Bflags1_1(:,f)=Bflags1;
               end
               
               AdflagH=sum(Ad_1(1:(ha-nha),:),1)/(ha-nha)*100;
               AdflagL=sum(Ad_1((ha-nha+1):end,:),1)/(la-nla)*100;
               
               Pflag1H=sum(Pflags1_1(1:(gp-ngp),:),1)/(gp-ngp)*100;
               Aflag1H=sum(Aflags1_1(1:(ga-nga),:),1)/(ga-nga)*100;
               Cflag1H=sum(Cflags1_1(1:(gc-ngc),:),1)/(gc-ngc)*100;
               Eflag1H=sum(Eflags1_1(1:(ge-nge),:),1)/(ge-nge)*100;
               Bflag1H=sum(Bflags1_1(1:(gb-ngb),:),1)/(gb-ngb)*100;
               
               Pflag1L=sum(Pflags1_1((gp-ngp+1):end,:),1)/(bp-nbp)*100;
               Aflag1L=sum(Aflags1_1((ga-nga+1):end,:),1)/(ba-nba)*100;
               Cflag1L=sum(Cflags1_1((gc-ngc+1):end,:),1)/(bc-nbc)*100;
               Eflag1L=sum(Eflags1_1((ge-nge+1):end,:),1)/(be-nbe)*100;
               Bflag1L=sum(Bflags1_1((gb-ngb+1):end,:),1)/(bb-nbb)*100;
               
               avergha=sum(Adflags(1:ha,:),1)/ha*100;
               avergla=sum(Adflags(ha+1:end,:),1)/la*100;
               averghcp=sum(Pflags(1:gp,:),1)/(gp)*100;
               averghca=sum(Aflags(1:ga,:),1)/(ga)*100;
               averghcc=sum(Cflags(1:gc,:),1)/(gc)*100;
               averghce=sum(Eflags(1:ge,:),1)/(ge)*100;
               averghcb=sum(Bflags(1:gb,:),1)/(gb)*100;
               averglcp=sum(Pflags(gp+1:end,:),1)/(bp)*100;
               averglca=sum(Aflags(ga+1:end,:),1)/(ba)*100;
               averglcc=sum(Cflags(gc+1:end,:),1)/(bc)*100;
               averglce=sum(Eflags(ge+1:end,:),1)/(be)*100;
               averglcb=sum(Bflags(gb+1:end,:),1)/(bb)*100;
              
               stdha=std(AdflagH);
               stdla=std(AdflagL);
               stdhcp=std(Pflag1H);
               stdlcp=std(Pflag1L);
               stdhca=std(Aflag1H);
               stdlca=std(Aflag1L);
               stdhcc=std(Cflag1H);
               stdlcc=std(Cflag1L);
               stdhce=std(Eflag1H);
               stdlce=std(Eflag1L);
               stdhcb=std(Bflag1H);
               stdlcb=std(Bflag1L);  
               if pvalue<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10 
                  opt{x+13}(q2,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
                  opt_0=zeros(1,13-x);
                  opt_1{x+13}(q2,:)=[100,opt_0,Polar_num(Polar1{x}(i,:))];
                  opt_2{x+13}(q2,:)=Flags;
                  q2=q2+1;
               elseif pvalue1<0.05 &&averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10  
                  opt{x+13}(q2,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
                  opt_0=zeros(1,13-x);
                  opt_1{x+13}(q2,:)=[100,opt_0,Polar_num(Polar1{x}(i,:))];
                  opt_2{x+13}(q2,:)=Flags;
                  q2=q2+1;
               end
            end
        end
    end
end
%% single amino acid
Single=netcharge;
PSingle=Pnetcharge;
ASingle=Anetcharge;
CSingle=Cnetcharge;
ESingle=Enetcharge;
BSingle=Bnetcharge;
AdSingle=ABnetcharge;
y=1;
for a=0:1:max_nc
    flags=Single>a;
    for f=1:50
        Fl=flags*1;
        Fl(clim(:,f),:)=[];
        Fl_1(:,f)=Fl;
    end
    averghc=(sum(flags(1:(hc),:),1)/(hc)*100);
            averglc=(sum(flags((hc+1):end,:),1)/(lc)*100);
            stdhc=std(sum(Fl_1(1:(hc-nhc)),1)/(hc-nhc)*100);
            stdlc=std(sum(Fl_1((hc-nhc+1):end),1)/(lc-nlc)*100);
            if averghc-stdhc<high && averglc+stdlc>low
               mAb=sum(flags)/(hc+lc)*100;
               highspec=averghc;
               lowspec=averglc;
               flag=a;
               Flags=flags';
               table(1,1)=sum(flags(1:hc)<1);
               table(2,1)=hc-table(1,1);
               table(1,2)=sum(flags(hc+1:end)<1);
               table(2,2)=lc-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=hc+lc;
               pvalue=TwoTailed(table);
       
       Adflags=AdSingle>a;
       Pflags=PSingle>a;
       Aflags=ASingle>a;
       Cflags=CSingle>a;
       Eflags=ESingle>a;
       Bflags=BSingle>a;
       
       table(1,1)=sum(Adflags(1:ha)<1);
               table(2,1)=ha-table(1,1);
               table(1,2)=sum(Adflags(ha+1:end)<1);
               table(2,2)=la-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=ha+la;
               pvalue1=TwoTailed(table);
       for f=1:50
                   Ad=Adflags;
                   Pflags1=Pflags;
                   Aflags1=Aflags;
                   Cflags1=Cflags;
                   Eflags1=Eflags;
                   Bflags1=Bflags;
                   
                   Ad(adi(:,f),:)=[];
                   Pflags1(clis_p(:,f),:)=[];
                   Aflags1(clis_a(:,f),:)=[];
                   Cflags1(clis_c(:,f),:)=[];
                   Eflags1(clis_e(:,f),:)=[];
                   Bflags1(clis_b(:,f),:)=[];
                   
                   Ad_1(:,f)=Ad;
                   Pflags1_1(:,f)=Pflags1;
                   Aflags1_1(:,f)=Aflags1;
                   Cflags1_1(:,f)=Cflags1;
                   Eflags1_1(:,f)=Eflags1;
                   Bflags1_1(:,f)=Bflags1;
               end
               
               AdflagH=sum(Ad_1(1:(ha-nha),:),1)/(ha-nha)*100;
               AdflagL=sum(Ad_1((ha-nha+1):end,:),1)/(la-nla)*100;
               
               Pflag1H=sum(Pflags1_1(1:(gp-ngp),:),1)/(gp-ngp)*100;
               Aflag1H=sum(Aflags1_1(1:(ga-nga),:),1)/(ga-nga)*100;
               Cflag1H=sum(Cflags1_1(1:(gc-ngc),:),1)/(gc-ngc)*100;
               Eflag1H=sum(Eflags1_1(1:(ge-nge),:),1)/(ge-nge)*100;
               Bflag1H=sum(Bflags1_1(1:(gb-ngb),:),1)/(gb-ngb)*100;
               
               Pflag1L=sum(Pflags1_1((gp-ngp+1):end,:),1)/(bp-nbp)*100;
               Aflag1L=sum(Aflags1_1((ga-nga+1):end,:),1)/(ba-nba)*100;
               Cflag1L=sum(Cflags1_1((gc-ngc+1):end,:),1)/(bc-nbc)*100;
               Eflag1L=sum(Eflags1_1((ge-nge+1):end,:),1)/(be-nbe)*100;
               Bflag1L=sum(Bflags1_1((gb-ngb+1):end,:),1)/(bb-nbb)*100;
               
              avergha=sum(Adflags(1:ha,:),1)/ha*100;
               avergla=sum(Adflags(ha+1:end,:),1)/la*100;
               averghcp=sum(Pflags(1:gp,:),1)/(gp)*100;
               averghca=sum(Aflags(1:ga,:),1)/(ga)*100;
               averghcc=sum(Cflags(1:gc,:),1)/(gc)*100;
               averghce=sum(Eflags(1:ge,:),1)/(ge)*100;
               averghcb=sum(Bflags(1:gb,:),1)/(gb)*100;
               averglcp=sum(Pflags(gp+1:end,:),1)/(bp)*100;
               averglca=sum(Aflags(ga+1:end,:),1)/(ba)*100;
               averglcc=sum(Cflags(gc+1:end,:),1)/(bc)*100;
               averglce=sum(Eflags(ge+1:end,:),1)/(be)*100;
               averglcb=sum(Bflags(gb+1:end,:),1)/(bb)*100;
              
               stdha=std(AdflagH);
               stdla=std(AdflagL);
               stdhcp=std(Pflag1H);
               stdlcp=std(Pflag1L);
               stdhca=std(Aflag1H);
               stdlca=std(Aflag1L);
               stdhcc=std(Cflag1H);
               stdlcc=std(Cflag1L);
               stdhce=std(Eflag1H);
               stdlce=std(Eflag1L);
               stdhcb=std(Bflag1H);
               stdlcb=std(Bflag1L);
       if pvalue<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10  
          opt{45}(y,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
          opt_1{45}(y,:)=[710,0,0,0,0,0,0,0,0,0,0,0,0,1];
          opt_2{45}(y,:)=[Flags];
          y=y+1;
       elseif pvalue1<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10  
          opt{45}(y,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
          opt_1{45}(y,:)=[710,0,0,0,0,0,0,0,0,0,0,0,0,1];
          opt_2{45}(y,:)=[Flags];
          y=y+1;
       end
    end
end
               
exist opt               
if ans~=0
   [~,r]=size(opt);
   Res=[];
   for i=1:r
       Res=[Res;opt_1{i} opt{i} opt_2{i}];
   end
   Res(isnan(Res(:,1)),:) = [] ;

   [e,~]=size(Res);
   Fx=Res(:,1:14);
   Px=Res(:,15:19);
   q=0;
   for i=1:e-1
       for j=i+1:e
           p=ismember(Fx(i,:),Fx(j,:));
           p1=sum(p);
           if p1==14 && Px(i,1)==Px(j,1)&& Px(i,2)==Px(j,2) && Px(i,3)==Px(j,3)
              q=q+1;
              Row(q)=j;
           end
       end
   end

   if q~=0
      Row=unique(Row);
      Res(Row,:) = []; 
      Res(:,18)=-1*Res(:,18);
      Res(:,20)=-1*Res(:,20);
      Res(:,22)=-1*Res(:,22);
      Res(:,24)=-1*Res(:,24);
      Res(:,26)=-1*Res(:,26);
      Res(:,28)=-1*Res(:,28);
      Res(:,30)=-1*Res(:,30);
      Res=sortrows(Res,[16,18,30,20,22,24,26,28]);
      Res(:,18)=-1*Res(:,18);
      Res(:,20)=-1*Res(:,20);
      Res(:,22)=-1*Res(:,22);
      Res(:,24)=-1*Res(:,24);
      Res(:,26)=-1*Res(:,26);
      Res(:,28)=-1*Res(:,28);
      Res(:,30)=-1*Res(:,30);
      [res,~]=size(Res);
      Sig=zeros(res,1);
      for i=1:res
          pc=i/res*0.05;
          if Res(i,16)<pc
             Sig(i,:)=1;
          end
      end
      s=find(Sig==1);
      [rs,~]=size(s);
      if rs==0
         disp('No results')
      else
      max1=max(s);
      Res_1=Res(1:max1,:); % results from FDR
      end
   else
      Res(:,18)=-1*Res(:,18);
      Res(:,20)=-1*Res(:,20);
      Res(:,22)=-1*Res(:,22);
      Res(:,24)=-1*Res(:,24);
      Res(:,26)=-1*Res(:,26);
      Res(:,28)=-1*Res(:,28);
      Res(:,30)=-1*Res(:,30);
      Res=sortrows(Res,[16,18,30,20,22,24,26,28]);
      Res(:,18)=-1*Res(:,18);
      Res(:,20)=-1*Res(:,20);
      Res(:,22)=-1*Res(:,22);
      Res(:,24)=-1*Res(:,24);
      Res(:,26)=-1*Res(:,26);
      Res(:,28)=-1*Res(:,28);
      Res(:,30)=-1*Res(:,30);
      [res,~]=size(Res);
      Sig=zeros(res,1);
      for i=1:res
          pc=i/res*0.05;
          if Res(i,16)<pc
          Sig(i,:)=1;
          end
      end
      s=find(Sig==1);
      [rs,~]=size(s);
      if rs==0
         disp('No results')
      else
         max1=max(s);
         Res_1=Res(1:max1,:); % results from FDR
      end
   end
   exist Res_1
   if ans ~=0 
      [res1,~]=size(Res_1);
      if res1~=0
         Sign_v=Res_1(:,1);
         Flag_AA=Res_1(:,2:14);
         other_info=Res_1(:,15:end);
         for i=1:res1
             F=nonzeros(Flag_AA(i,:))';
             if Sign_v(i)==510
                Hydro_AA='GAVLIMFWPYRKH';
                Flag_AA1{i}=Hydro_AA(F);
                Sign1{i}='>';
             elseif Sign_v(i)==100
                Pola_AA='GASTQNYDE';
                Flag_AA1{i}=Pola_AA(F);
                Sign1{i}='<';
             elseif Sign_v(i)==710
                Flag_AA1{i}='charge';
                Sign1{i}='>';
             end
         end
      end
      Flag=cell2table(Flag_AA1');
      Sign=cell2table(Sign1');
      Flag.Properties.VariableNames{'Var1'}='Flag';
      Sign.Properties.VariableNames{'Var1'}='Sign';
      info=array2table(other_info);
      R=cell(res1,1);
      R(:)={Region(re2)};
      Singleflag=[R,Flag,Sign,info];
   else
      Singleflag=[];
   end
   
  else
   Singleflag=[];
end
re2
  Finalflag=[Finalflag;Singleflag];
end
Region1={'H1','H1','H1','H1','H1','H1','H1','H1','H1','H1','H2','H2','H2','H2','H2','H2','H3','H3','H3','L1'};
Region2={'H2','H2','H2','H2','H3','H3','H3','L1','L1','L2','H3','H3','H3','L1','L1','L2','L1','L1','L2','L2'};
Region3={'H3','L1','L2','L3','L1','L2','L3','L2','L3','L3','L1','L2','L3','L2','L3','L3','L2','L3','L3','L3'};
Region={'H1H2H3','H1H2L1','H1H2L2','H1H2L3','H1H3L1','H1H3L2','H1H3L3','H1L1L2','H1L1L3','H1L2L3','H2H3L1','H2H3L2','H2H3L3','H2L1L2','H2L1L3','H2L2L3','H3L1L2','H3L1L3','H3L2L3','L1L2L3'};
for re3=1:20
    clear opt opt_0 opt_1 opt_2 Flag_AA1 Sign1 Sign_v Res_1 Row
warning('off','all');
CAB1=xlsread('Clinical_AA',Region1{re3},'A1:T200');
%HAB=xlsread('human antibody_sequence_AA comp_newlimit_v1',Region,'A2:T20013');
CP1=xlsread('Clinical_AA_PSR',Region1{re3},'A1:T300');
CA1=xlsread('Clinical_AA_ACSINS',Region1{re3},'A1:T300');
CC1=xlsread('Clinical_AA_CSI',Region1{re3},'A1:T300');
CE1=xlsread('Clinical_AA_ELISA',Region1{re3},'A1:T300');
CB1=xlsread('Clinical_AA_BVP',Region1{re3},'A1:T300');
AAB1=xlsread('Human_AA_100sim',Region1{re3},'A1:T2000');

CAB2=xlsread('Clinical_AA',Region2{re3},'A1:T200');
%HAB=xlsread('human antibody_sequence_AA comp_newlimit_v1',Region,'A2:T20013');
CP2=xlsread('Clinical_AA_PSR',Region2{re3},'A1:T300');
CA2=xlsread('Clinical_AA_ACSINS',Region2{re3},'A1:T300');
CC2=xlsread('Clinical_AA_CSI',Region2{re3},'A1:T300');
CE2=xlsread('Clinical_AA_ELISA',Region2{re3},'A1:T300');
CB2=xlsread('Clinical_AA_BVP',Region2{re3},'A1:T300');
AAB2=xlsread('Human_AA_100sim',Region2{re3},'A1:T2000');


CAB3=xlsread('Clinical_AA',Region3{re3},'A1:T200');
%HAB=xlsread('human antibody_sequence_AA comp_newlimit_v1',Region,'A2:T20013');
CP3=xlsread('Clinical_AA_PSR',Region3{re3},'A1:T300');
CA3=xlsread('Clinical_AA_ACSINS',Region3{re3},'A1:T300');
CC3=xlsread('Clinical_AA_CSI',Region3{re3},'A1:T300');
CE3=xlsread('Clinical_AA_ELISA',Region3{re3},'A1:T300');
CB3=xlsread('Clinical_AA_BVP',Region3{re3},'A1:T300');
AAB3=xlsread('Human_AA_100sim',Region3{re3},'A1:T2000');

CAB=CAB1+CAB2+CAB3;
CP=CP1+CP2+CP3;
CA=CA1+CA2+CA3;
CC=CC1+CC2+CC3;
CE=CE1+CE2+CE3;
CB=CB1+CB2+CB3;
AAB=AAB1+AAB2+AAB3;
adihs=[];adils=[];clis_ph=[];clis_pl=[];clis_ch=[];clis_cl=[];clis_eh=[];clis_el=[];clis_bh=[];clis_bl=[];
clis_ah=[];clis_al=[];clihs=[];clils=[];

adi=xlsread('2ndset5_combine_singleflags_v2_new','adimab','A2:CV100');
clis_p=xlsread('2ndset5_combine_singleflags_v2_new','PSR','A2:CV20');
clis_a=xlsread('2ndset5_combine_singleflags_v2_new','ACSINS','A2:CV20');
clis_c=xlsread('2ndset5_combine_singleflags_v2_new','CSI','A2:CV20');
clis_e=xlsread('2ndset5_combine_singleflags_v2_new','ELISA','A2:CV20');
clis_b=xlsread('2ndset5_combine_singleflags_v2_new','BVP','A2:CV20');
clim=xlsread('2ndset5_combine_singleflags_v2_new','clinical','A2:CV100');
ngp=sum(clis_p(:,1)<=gp);nbp=sum(clis_p(:,1)>gp);
nga=sum(clis_a(:,1)<=ga);nba=sum(clis_a(:,1)>ga);
ngc=sum(clis_c(:,1)<=gc);nbc=sum(clis_c(:,1)>gc);
nge=sum(clis_e(:,1)<=ge);nbe=sum(clis_e(:,1)>ge);
ngb=sum(clis_b(:,1)<=gb);nbb=sum(clis_b(:,1)>gb);
nha=sum(adi(:,1)<=ha);nla=sum(adi(:,1)>ha);
nhc=sum(clim(:,1)<=hc);nlc=sum(clim(:,1)>hc);
%hv=1033;lv=59;hr=104;lr=6;highP=0.15;
A=CAB(Total2,1);C=CAB(Total2,2);D=CAB(Total2,3);E=CAB(Total2,4);F=CAB(Total2,5);G=CAB(Total2,6);H=CAB(Total2,7);
I=CAB(Total2,8);K=CAB(Total2,9);L=CAB(Total2,10);M=CAB(Total2,11);N=CAB(Total2,12);P=CAB(Total2,13);Q=CAB(Total2,14);
R=CAB(Total2,15);S=CAB(Total2,16);T=CAB(Total2,17);V=CAB(Total2,18);W=CAB(Total2,19);Y=CAB(Total2,20);
netcharge=R+K+0.1*H-D-E;
max_nc=floor(max(netcharge));
%% Human Antibody information PSR
A1=CP(PSR2,1);C1=CP(PSR2,2);D1=CP(PSR2,3);E1=CP(PSR2,4);F1=CP(PSR2,5);G1=CP(PSR2,6);H1=CP(PSR2,7);
I1=CP(PSR2,8);K1=CP(PSR2,9);L1=CP(PSR2,10);M1=CP(PSR2,11);N1=CP(PSR2,12);P1=CP(PSR2,13);Q1=CP(PSR2,14);
R1=CP(PSR2,15);S1=CP(PSR2,16);T1=CP(PSR2,17);V1=CP(PSR2,18);W1=CP(PSR2,19);Y1=CP(PSR2,20);
Pnetcharge=R1+K1+0.1*H1-D1-E1;
%% Human Antibody information ACSINS
A2=CA(ACSINS2,1);C2=CA(ACSINS2,2);D2=CA(ACSINS2,3);E2=CA(ACSINS2,4);F2=CA(ACSINS2,5);G2=CA(ACSINS2,6);H2=CA(ACSINS2,7);
I2=CA(ACSINS2,8);K2=CA(ACSINS2,9);L2=CA(ACSINS2,10);M2=CA(ACSINS2,11);N2=CA(ACSINS2,12);P2=CA(ACSINS2,13);Q2=CA(ACSINS2,14);
R2=CA(ACSINS2,15);S2=CA(ACSINS2,16);T2=CA(ACSINS2,17);V2=CA(ACSINS2,18);W2=CA(ACSINS2,19);Y2=CA(ACSINS2,20);
Anetcharge=R2+K2+0.1*H2-D2-E2;
%% Human Antibody information CSI
A3=CC(CSI2,1);D3=CC(CSI2,3);E3=CC(CSI2,4);F3=CC(CSI2,5);G3=CC(CSI2,6);H3=CC(CSI2,7);
I3=CC(CSI2,8);K3=CC(CSI2,9);L3=CC(CSI2,10);M3=CC(CSI2,11);N3=CC(CSI2,12);P3=CC(CSI2,13);Q3=CC(CSI2,14);
R3=CC(CSI2,15);S3=CC(CSI2,16);T3=CC(CSI2,17);V3=CC(CSI2,18);W3=CC(CSI2,19);Y3=CC(CSI2,20);
Cnetcharge=R3+K3+0.1*H3-D3-E3;
%% Human Antibody information ELISA
A4=CE(ELISA2,1);C4=CE(ELISA2,2);D4=CE(ELISA2,3);E4=CE(ELISA2,4);F4=CE(ELISA2,5);G4=CE(ELISA2,6);H4=CE(ELISA2,7);
I4=CE(ELISA2,8);K4=CE(ELISA2,9);L4=CE(ELISA2,10);M4=CE(ELISA2,11);N4=CE(ELISA2,12);P4=CE(ELISA2,13);Q4=CE(ELISA2,14);
R4=CE(ELISA2,15);S4=CE(ELISA2,16);T4=CE(ELISA2,17);V4=CE(ELISA2,18);W4=CE(ELISA2,19);Y4=CE(ELISA2,20);
Enetcharge=R4+K4+0.1*H4-D4-E4;
%% Human Antibody information ELISA
A5=CB(BVP2,1);C5=CB(BVP2,2);D5=CB(BVP2,3);E5=CB(BVP2,4);F5=CB(BVP2,5);G5=CB(BVP2,6);H5=CB(BVP2,7);
I5=CB(BVP2,8);K5=CB(BVP2,9);L5=CB(BVP2,10);M5=CB(BVP2,11);N5=CB(BVP2,12);P5=CB(BVP2,13);Q5=CB(BVP2,14);
R5=CB(BVP2,15);S5=CB(BVP2,16);T5=CB(BVP2,17);V5=CB(BVP2,18);W5=CB(BVP2,19);Y5=CB(BVP2,20);
Bnetcharge=R5+K5+0.1*H5-D5-E5;
%% Human B cell from Adimab
ABA=AAB(train2,1);ABC=AAB(train2,2);ABD=AAB(train2,3);ABE=AAB(train2,4);ABF=AAB(train2,5);ABG=AAB(train2,6);ABH=AAB(train2,7);
ABI=AAB(train2,8);ABK=AAB(train2,9);ABL=AAB(train2,10);ABM=AAB(train2,11);ABN=AAB(train2,12);ABP=AAB(train2,13);ABQ=AAB(train2,14);
ABR=AAB(train2,15);ABS=AAB(train2,16);ABT=AAB(train2,17);ABV=AAB(train2,18);ABW=AAB(train2,19);ABY=AAB(train2,20);
ABnetcharge=ABR+ABK+0.1*ABH-ABD-ABE;

%% combination
% hydrophybic residue GAVLIMFWPYRKH
Hydro=[G,A,V,L,I,M,F,W,P,Y,R,K,H];
PHydro=[G1,A1,V1,L1,I1,M1,F1,W1,P1,Y1,R1,K1,H1];
AHydro=[G2,A2,V2,L2,I2,M2,F2,W2,P2,Y2,R2,K2,H2];
CHydro=[G3,A3,V3,L3,I3,M3,F3,W3,P3,Y3,R3,K3,H3];
EHydro=[G4,A4,V4,L4,I4,M4,F4,W4,P4,Y4,R4,K4,H4];
BHydro=[G5,A5,V5,L5,I5,M5,F5,W5,P5,Y5,R5,K5,H5];
AdHydro=[ABG,ABA,ABV,ABL,ABI,ABM,ABF,ABW,ABP,ABY,ABR,ABK,ABH];
Hydro_num=[1,2,3,4,5,6,7,8,9,10,11,12,13];
for i=1:13
    Hydro1{i}=combntns(1:13,i);
    [r(i),~]=size(Hydro1{i});
end    
for k=1:13
    for i=1:1:r(k)
        Hydro_1{k}(:,i)=sum(Hydro(:,Hydro1{k}(i,:)),2);
        PHydro_1{k}(:,i)=sum(PHydro(:,Hydro1{k}(i,:)),2);
        AHydro_1{k}(:,i)=sum(AHydro(:,Hydro1{k}(i,:)),2);
        CHydro_1{k}(:,i)=sum(CHydro(:,Hydro1{k}(i,:)),2);
        EHydro_1{k}(:,i)=sum(EHydro(:,Hydro1{k}(i,:)),2);
        BHydro_1{k}(:,i)=sum(BHydro(:,Hydro1{k}(i,:)),2);
        AdHydro_1{k}(:,i)=sum(AdHydro(:,Hydro1{k}(i,:)),2);
    end
    [~,r_1(k)]=size(Hydro_1{k});
    for j=1:1:r_1(k)
        maxbind_hy{k}(j)=max(Hydro_1{k}(:,j));
        minbind_hy{k}(j)=min(Hydro_1{k}(:,j));
        hy{k}(j)=(maxbind_hy{k}(j)-minbind_hy{k}(j))+1;
    end        
end    
for x=1:1:13
    q1=1;
    opt_0=[];
    for i=1:1:r(x)
        for a=minbind_hy{x}(i):1:maxbind_hy{x}(i)
            flags=Hydro_1{x}(:,i)>a;
            for f=1:50
                Fl=flags*1;
                Fl(clim(:,f),:)=[];
                Fl_1(:,f)=Fl;
            end
            averghc=(sum(flags(1:(hc),:),1)/(hc)*100);
            averglc=(sum(flags((hc+1):end,:),1)/(lc)*100);
            stdhc=std(sum(Fl_1(1:(hc-nhc)),1)/(hc-nhc)*100);
            stdlc=std(sum(Fl_1((hc-nhc+1):end),1)/(lc-nlc)*100);
            if averghc-stdhc<high && averglc+stdlc>low
               mAb=sum(flags)/(hc+lc)*100;
               highspec=averghc;
               lowspec=averglc;
               flag=a;
               Flags=flags';
               table(1,1)=sum(flags(1:hc)<1);
               table(2,1)=hc-table(1,1);
               table(1,2)=sum(flags(hc+1:end)<1);
               table(2,2)=lc-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=hc+lc;
               pvalue=TwoTailed(table);
               
               Adflags=AdHydro_1{x}(:,i)>a;
               Pflags=PHydro_1{x}(:,i)>a;
               Aflags=AHydro_1{x}(:,i)>a;
               Cflags=CHydro_1{x}(:,i)>a;
               Eflags=EHydro_1{x}(:,i)>a;
               Bflags=BHydro_1{x}(:,i)>a;
               
               table(1,1)=sum(Adflags(1:ha)<1);
               table(2,1)=ha-table(1,1);
               table(1,2)=sum(Adflags(ha+1:end)<1);
               table(2,2)=la-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=ha+la;
               pvalue1=TwoTailed(table);
               
               
               for f=1:50
                   Ad=Adflags;
                   Pflags1=Pflags;
                   Aflags1=Aflags;
                   Cflags1=Cflags;
                   Eflags1=Eflags;
                   Bflags1=Bflags;
                   
                   Ad(adi(:,f),:)=[];
                   Pflags1(clis_p(:,f),:)=[];
                   Aflags1(clis_a(:,f),:)=[];
                   Cflags1(clis_c(:,f),:)=[];
                   Eflags1(clis_e(:,f),:)=[];
                   Bflags1(clis_b(:,f),:)=[];
                   
                   Ad_1(:,f)=Ad;
                   Pflags1_1(:,f)=Pflags1;
                   Aflags1_1(:,f)=Aflags1;
                   Cflags1_1(:,f)=Cflags1;
                   Eflags1_1(:,f)=Eflags1;
                   Bflags1_1(:,f)=Bflags1;
               end
               
               AdflagH=sum(Ad_1(1:(ha-nha),:),1)/(ha-nha)*100;
               AdflagL=sum(Ad_1((ha-nha+1):end,:),1)/(la-nla)*100;
               
               Pflag1H=sum(Pflags1_1(1:(gp-ngp),:),1)/(gp-ngp)*100;
               Aflag1H=sum(Aflags1_1(1:(ga-nga),:),1)/(ga-nga)*100;
               Cflag1H=sum(Cflags1_1(1:(gc-ngc),:),1)/(gc-ngc)*100;
               Eflag1H=sum(Eflags1_1(1:(ge-nge),:),1)/(ge-nge)*100;
               Bflag1H=sum(Bflags1_1(1:(gb-ngb),:),1)/(gb-ngb)*100;
               
               Pflag1L=sum(Pflags1_1((gp-ngp+1):end,:),1)/(bp-nbp)*100;
               Aflag1L=sum(Aflags1_1((ga-nga+1):end,:),1)/(ba-nba)*100;
               Cflag1L=sum(Cflags1_1((gc-ngc+1):end,:),1)/(bc-nbc)*100;
               Eflag1L=sum(Eflags1_1((ge-nge+1):end,:),1)/(be-nbe)*100;
               Bflag1L=sum(Bflags1_1((gb-ngb+1):end,:),1)/(bb-nbb)*100;
               
               avergha=sum(Adflags(1:ha,:),1)/ha*100;
               avergla=sum(Adflags(ha+1:end,:),1)/la*100;
               averghcp=sum(Pflags(1:gp,:),1)/(gp)*100;
               averghca=sum(Aflags(1:ga,:),1)/(ga)*100;
               averghcc=sum(Cflags(1:gc,:),1)/(gc)*100;
               averghce=sum(Eflags(1:ge,:),1)/(ge)*100;
               averghcb=sum(Bflags(1:gb,:),1)/(gb)*100;
               averglcp=sum(Pflags(gp+1:end,:),1)/(bp)*100;
               averglca=sum(Aflags(ga+1:end,:),1)/(ba)*100;
               averglcc=sum(Cflags(gc+1:end,:),1)/(bc)*100;
               averglce=sum(Eflags(ge+1:end,:),1)/(be)*100;
               averglcb=sum(Bflags(gb+1:end,:),1)/(bb)*100;
               stdha=std(AdflagH);
               stdla=std(AdflagL);
               stdhcp=std(Pflag1H);
               stdlcp=std(Pflag1L);
               stdhca=std(Aflag1H);
               stdlca=std(Aflag1L);
               stdhcc=std(Cflag1H);
               stdlcc=std(Cflag1L);
               stdhce=std(Eflag1H);
               stdlce=std(Eflag1L);
               stdhcb=std(Bflag1H);
               stdlcb=std(Bflag1L);
               
                  
               if pvalue<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10 
                  opt{x}(q1,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
                  opt_0=zeros(1,13-x);
                  opt_1{x}(q1,:)=[510,opt_0,Hydro_num(Hydro1{x}(i,:))];
                  opt_2{x}(q1,:)=Flags;
                  q1=q1+1;
                  
               elseif pvalue1<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10  
                  opt{x}(q1,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
                  opt_0=zeros(1,13-x);
                  opt_1{x}(q1,:)=[510,opt_0,Hydro_num(Hydro1{x}(i,:))];
                  opt_2{x}(q1,:)=Flags;
                  q1=q1+1;
               end
            end
        end
    end
end
%% combination 
% polar residue GAVLIMFWP
Polar=[G,A,S,T,Q,N,Y,D,E];
PPolar=[G1,A1,S1,T1,Q1,N1,Y1,D1,E1];
APolar=[G2,A2,S2,T2,Q2,N2,Y2,D2,E2];
CPolar=[G3,A3,S3,T3,Q3,N3,Y3,D3,E3];
EPolar=[G4,A4,S4,T4,Q4,N4,Y4,D4,E4];
BPolar=[G5,A5,S5,T5,Q5,N5,Y5,D5,E5];
AdPolar=[ABG,ABA,ABS,ABT,ABQ,ABN,ABY,ABD,ABE];
Polar_num=[1,2,3,4,5,6,7,8,9];
for i=1:9
    Polar1{i}=combntns(1:9,i);
    [r_p(i),~]=size(Polar1{i});
end 
for k=1:9
    for i=1:1:r_p(k)
        Polar_1{k}(:,i)=sum(Polar(:,Polar1{k}(i,:)),2);
        PPolar_1{k}(:,i)=sum(PPolar(:,Polar1{k}(i,:)),2);
        APolar_1{k}(:,i)=sum(APolar(:,Polar1{k}(i,:)),2);
        CPolar_1{k}(:,i)=sum(CPolar(:,Polar1{k}(i,:)),2);
        EPolar_1{k}(:,i)=sum(EPolar(:,Polar1{k}(i,:)),2);
        BPolar_1{k}(:,i)=sum(BPolar(:,Polar1{k}(i,:)),2);
        AdPolar_1{k}(:,i)=sum(AdPolar(:,Polar1{k}(i,:)),2);
    end
    [~,r_p1(k)]=size(Polar_1{k});
    for j=1:1:r_p1(k)
        maxbind_po{k}(j)=max(Polar_1{k}(:,j));
        minbind_po{k}(j)=min(Polar_1{k}(:,j));
        po{k}(j)=(maxbind_po{k}(j)-minbind_po{k}(j))+1;
    end        
end

for x=1:1:9
    q2=1;
    opt_0=[];
    for i=1:1:r_p(x)
        for a=minbind_po{x}(i):1:maxbind_po{x}(i)
             flags=Polar_1{x}(:,i)<a;
            for f=1:50
                Fl=flags*1;
                Fl(clim(:,f),:)=[];
                Fl_1(:,f)=Fl;
            end
            averghc=sum(flags(1:(hc),:),1)/(hc)*100;
            averglc=(sum(flags((hc+1):end,:),1)/(lc)*100);
            stdhc=std(sum(Fl_1(1:(hc-nhc)),1)/(hc-nhc)*100);
            stdlc=std(sum(Fl_1((hc-nhc+1):end),1)/(lc-nlc)*100);
            if averghc-stdhc<high && averglc+stdlc>low
               mAb=sum(flags)/(hc+lc)*100;
               highspec=averghc;
               lowspec=averglc;
               flag=a;
               Flags=flags';
               table(1,1)=sum(flags(1:hc)<1);
               table(2,1)=hc-table(1,1);
               table(1,2)=sum(flags(hc+1:end)<1);
               table(2,2)=lc-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=hc+lc;
               pvalue=TwoTailed(table);
               
               Adflags=AdPolar_1{x}(:,i)<a;
               Pflags=PPolar_1{x}(:,i)<a;
               Aflags=APolar_1{x}(:,i)<a;
               Cflags=CPolar_1{x}(:,i)<a;
               Eflags=EPolar_1{x}(:,i)<a;
               Bflags=BPolar_1{x}(:,i)<a;
               
               table(1,1)=sum(Adflags(1:ha)<1);
               table(2,1)=ha-table(1,1);
               table(1,2)=sum(Adflags(ha+1:end)<1);
               table(2,2)=la-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=ha+la;
               pvalue1=TwoTailed(table);
               
               for f=1:50
                   Ad=Adflags;
                   Pflags1=Pflags;
                   Aflags1=Aflags;
                   Cflags1=Cflags;
                   Eflags1=Eflags;
                   Bflags1=Bflags;
                   
                   Ad(adi(:,f),:)=[];
                   Pflags1(clis_p(:,f),:)=[];
                   Aflags1(clis_a(:,f),:)=[];
                   Cflags1(clis_c(:,f),:)=[];
                   Eflags1(clis_e(:,f),:)=[];
                   Bflags1(clis_b(:,f),:)=[];
                   
                   Ad_1(:,f)=Ad;
                   Pflags1_1(:,f)=Pflags1;
                   Aflags1_1(:,f)=Aflags1;
                   Cflags1_1(:,f)=Cflags1;
                   Eflags1_1(:,f)=Eflags1;
                   Bflags1_1(:,f)=Bflags1;
               end
               
               AdflagH=sum(Ad_1(1:(ha-nha),:),1)/(ha-nha)*100;
               AdflagL=sum(Ad_1((ha-nha+1):end,:),1)/(la-nla)*100;
               
               Pflag1H=sum(Pflags1_1(1:(gp-ngp),:),1)/(gp-ngp)*100;
               Aflag1H=sum(Aflags1_1(1:(ga-nga),:),1)/(ga-nga)*100;
               Cflag1H=sum(Cflags1_1(1:(gc-ngc),:),1)/(gc-ngc)*100;
               Eflag1H=sum(Eflags1_1(1:(ge-nge),:),1)/(ge-nge)*100;
               Bflag1H=sum(Bflags1_1(1:(gb-ngb),:),1)/(gb-ngb)*100;
               
               Pflag1L=sum(Pflags1_1((gp-ngp+1):end,:),1)/(bp-nbp)*100;
               Aflag1L=sum(Aflags1_1((ga-nga+1):end,:),1)/(ba-nba)*100;
               Cflag1L=sum(Cflags1_1((gc-ngc+1):end,:),1)/(bc-nbc)*100;
               Eflag1L=sum(Eflags1_1((ge-nge+1):end,:),1)/(be-nbe)*100;
               Bflag1L=sum(Bflags1_1((gb-ngb+1):end,:),1)/(bb-nbb)*100;
               
               avergha=sum(Adflags(1:ha,:),1)/ha*100;
               avergla=sum(Adflags(ha+1:end,:),1)/la*100;
               averghcp=sum(Pflags(1:gp,:),1)/(gp)*100;
               averghca=sum(Aflags(1:ga,:),1)/(ga)*100;
               averghcc=sum(Cflags(1:gc,:),1)/(gc)*100;
               averghce=sum(Eflags(1:ge,:),1)/(ge)*100;
               averghcb=sum(Bflags(1:gb,:),1)/(gb)*100;
               averglcp=sum(Pflags(gp+1:end,:),1)/(bp)*100;
               averglca=sum(Aflags(ga+1:end,:),1)/(ba)*100;
               averglcc=sum(Cflags(gc+1:end,:),1)/(bc)*100;
               averglce=sum(Eflags(ge+1:end,:),1)/(be)*100;
               averglcb=sum(Bflags(gb+1:end,:),1)/(bb)*100;
              
               stdha=std(AdflagH);
               stdla=std(AdflagL);
               stdhcp=std(Pflag1H);
               stdlcp=std(Pflag1L);
               stdhca=std(Aflag1H);
               stdlca=std(Aflag1L);
               stdhcc=std(Cflag1H);
               stdlcc=std(Cflag1L);
               stdhce=std(Eflag1H);
               stdlce=std(Eflag1L);
               stdhcb=std(Bflag1H);
               stdlcb=std(Bflag1L);  
               if pvalue<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10 
                  opt{x+13}(q2,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
                  opt_0=zeros(1,13-x);
                  opt_1{x+13}(q2,:)=[100,opt_0,Polar_num(Polar1{x}(i,:))];
                  opt_2{x+13}(q2,:)=Flags;
                  q2=q2+1;
               elseif pvalue1<0.05 &&averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10  
                  opt{x+13}(q2,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
                  opt_0=zeros(1,13-x);
                  opt_1{x+13}(q2,:)=[100,opt_0,Polar_num(Polar1{x}(i,:))];
                  opt_2{x+13}(q2,:)=Flags;
                  q2=q2+1;
               end
            end
        end
    end
end

%% single amino acid
Single=netcharge;
PSingle=Pnetcharge;
ASingle=Anetcharge;
CSingle=Cnetcharge;
ESingle=Enetcharge;
BSingle=Bnetcharge;
AdSingle=ABnetcharge;
y=1;
for a=0:1:max_nc
    flags=Single>a;
    for f=1:50
        Fl=flags*1;
        Fl(clim(:,f),:)=[];
        Fl_1(:,f)=Fl;
    end
    averghc=(sum(flags(1:(hc),:),1)/(hc)*100);
            averglc=(sum(flags((hc+1):end,:),1)/(lc)*100);
            stdhc=std(sum(Fl_1(1:(hc-nhc)),1)/(hc-nhc)*100);
            stdlc=std(sum(Fl_1((hc-nhc+1):end),1)/(lc-nlc)*100);
            if averghc-stdhc<high && averglc+stdlc>low
               mAb=sum(flags)/(hc+lc)*100;
               highspec=averghc;
               lowspec=averglc;
               flag=a;
               Flags=flags';
               table(1,1)=sum(flags(1:hc)<1);
               table(2,1)=hc-table(1,1);
               table(1,2)=sum(flags(hc+1:end)<1);
               table(2,2)=lc-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=hc+lc;
               pvalue=TwoTailed(table);
       
       Adflags=AdSingle>a;
       Pflags=PSingle>a;
       Aflags=ASingle>a;
       Cflags=CSingle>a;
       Eflags=ESingle>a;
       Bflags=BSingle>a;
       
       table(1,1)=sum(Adflags(1:ha)<1);
               table(2,1)=ha-table(1,1);
               table(1,2)=sum(Adflags(ha+1:end)<1);
               table(2,2)=la-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=ha+la;
               pvalue1=TwoTailed(table);
       for f=1:50
                   Ad=Adflags;
                   Pflags1=Pflags;
                   Aflags1=Aflags;
                   Cflags1=Cflags;
                   Eflags1=Eflags;
                   Bflags1=Bflags;
                   
                   Ad(adi(:,f),:)=[];
                   Pflags1(clis_p(:,f),:)=[];
                   Aflags1(clis_a(:,f),:)=[];
                   Cflags1(clis_c(:,f),:)=[];
                   Eflags1(clis_e(:,f),:)=[];
                   Bflags1(clis_b(:,f),:)=[];
                   
                   Ad_1(:,f)=Ad;
                   Pflags1_1(:,f)=Pflags1;
                   Aflags1_1(:,f)=Aflags1;
                   Cflags1_1(:,f)=Cflags1;
                   Eflags1_1(:,f)=Eflags1;
                   Bflags1_1(:,f)=Bflags1;
               end
               
               AdflagH=sum(Ad_1(1:(ha-nha),:),1)/(ha-nha)*100;
               AdflagL=sum(Ad_1((ha-nha+1):end,:),1)/(la-nla)*100;
               
               Pflag1H=sum(Pflags1_1(1:(gp-ngp),:),1)/(gp-ngp)*100;
               Aflag1H=sum(Aflags1_1(1:(ga-nga),:),1)/(ga-nga)*100;
               Cflag1H=sum(Cflags1_1(1:(gc-ngc),:),1)/(gc-ngc)*100;
               Eflag1H=sum(Eflags1_1(1:(ge-nge),:),1)/(ge-nge)*100;
               Bflag1H=sum(Bflags1_1(1:(gb-ngb),:),1)/(gb-ngb)*100;
               
               Pflag1L=sum(Pflags1_1((gp-ngp+1):end,:),1)/(bp-nbp)*100;
               Aflag1L=sum(Aflags1_1((ga-nga+1):end,:),1)/(ba-nba)*100;
               Cflag1L=sum(Cflags1_1((gc-ngc+1):end,:),1)/(bc-nbc)*100;
               Eflag1L=sum(Eflags1_1((ge-nge+1):end,:),1)/(be-nbe)*100;
               Bflag1L=sum(Bflags1_1((gb-ngb+1):end,:),1)/(bb-nbb)*100;
               
               avergha=sum(Adflags(1:ha,:),1)/ha*100;
               avergla=sum(Adflags(ha+1:end,:),1)/la*100;
               averghcp=sum(Pflags(1:gp,:),1)/(gp)*100;
               averghca=sum(Aflags(1:ga,:),1)/(ga)*100;
               averghcc=sum(Cflags(1:gc,:),1)/(gc)*100;
               averghce=sum(Eflags(1:ge,:),1)/(ge)*100;
               averghcb=sum(Bflags(1:gb,:),1)/(gb)*100;
               averglcp=sum(Pflags(gp+1:end,:),1)/(bp)*100;
               averglca=sum(Aflags(ga+1:end,:),1)/(ba)*100;
               averglcc=sum(Cflags(gc+1:end,:),1)/(bc)*100;
               averglce=sum(Eflags(ge+1:end,:),1)/(be)*100;
               averglcb=sum(Bflags(gb+1:end,:),1)/(bb)*100;
              
               stdha=std(AdflagH);
               stdla=std(AdflagL);
               stdhcp=std(Pflag1H);
               stdlcp=std(Pflag1L);
               stdhca=std(Aflag1H);
               stdlca=std(Aflag1L);
               stdhcc=std(Cflag1H);
               stdlcc=std(Cflag1L);
               stdhce=std(Eflag1H);
               stdlce=std(Eflag1L);
               stdhcb=std(Bflag1H);
               stdlcb=std(Bflag1L);
       if pvalue<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10  
          opt{45}(y,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
          opt_1{45}(y,:)=[710,0,0,0,0,0,0,0,0,0,0,0,0,1];
          opt_2{45}(y,:)=[Flags];
          y=y+1;
       elseif pvalue1<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10  
          opt{45}(y,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
          opt_1{45}(y,:)=[710,0,0,0,0,0,0,0,0,0,0,0,0,1];
          opt_2{45}(y,:)=[Flags];
          y=y+1;
       end
    end
end
               
exist opt               
if ans~=0
   [~,r]=size(opt);
   Res=[];
   for i=1:r
       Res=[Res;opt_1{i} opt{i} opt_2{i}];
   end
   Res(isnan(Res(:,1)),:) = [] ;

   [e,~]=size(Res);
   Fx=Res(:,1:14);
   Px=Res(:,15:19);
   q=0;
   for i=1:e-1
       for j=i+1:e
           p=ismember(Fx(i,:),Fx(j,:));
           p1=sum(p);
           if p1==14 && Px(i,1)==Px(j,1)&& Px(i,2)==Px(j,2) && Px(i,3)==Px(j,3)
              q=q+1;
              Row(q)=j;
           end
       end
   end

   if q~=0
      Row=unique(Row);
      Res(Row,:) = []; 
      Res(:,18)=-1*Res(:,18);
      Res(:,20)=-1*Res(:,20);
      Res(:,22)=-1*Res(:,22);
      Res(:,24)=-1*Res(:,24);
      Res(:,26)=-1*Res(:,26);
      Res(:,28)=-1*Res(:,28);
      Res(:,30)=-1*Res(:,30);
      Res=sortrows(Res,[16,18,30,20,22,24,26,28]);
      Res(:,18)=-1*Res(:,18);
      Res(:,20)=-1*Res(:,20);
      Res(:,22)=-1*Res(:,22);
      Res(:,24)=-1*Res(:,24);
      Res(:,26)=-1*Res(:,26);
      Res(:,28)=-1*Res(:,28);
      Res(:,30)=-1*Res(:,30);
      [res,~]=size(Res);
      Sig=zeros(res,1);
      for i=1:res
          pc=i/res*0.05;
          if Res(i,16)<pc
             Sig(i,:)=1;
          end
      end
      s=find(Sig==1);
      [rs,~]=size(s);
      if rs==0
         disp('No results')
      else
      max1=max(s);
      Res_1=Res(1:max1,:); % results from FDR
      end
   else
      Res(:,18)=-1*Res(:,18);
      Res(:,20)=-1*Res(:,20);
      Res(:,22)=-1*Res(:,22);
      Res(:,24)=-1*Res(:,24);
      Res(:,26)=-1*Res(:,26);
      Res(:,28)=-1*Res(:,28);
      Res(:,30)=-1*Res(:,30);
      Res=sortrows(Res,[16,18,30,20,22,24,26,28]);
      Res(:,18)=-1*Res(:,18);
      Res(:,20)=-1*Res(:,20);
      Res(:,22)=-1*Res(:,22);
      Res(:,24)=-1*Res(:,24);
      Res(:,26)=-1*Res(:,26);
      Res(:,28)=-1*Res(:,28);
      Res(:,30)=-1*Res(:,30);
      [res,~]=size(Res);
      Sig=zeros(res,1);
      for i=1:res
          pc=i/res*0.05;
          if Res(i,16)<pc
          Sig(i,:)=1;
          end
      end
      s=find(Sig==1);
      [rs,~]=size(s);
      if rs==0
         disp('No results')
      else
         max1=max(s);
         Res_1=Res(1:max1,:); % results from FDR
      end
   end
   exist Res_1
   if ans ~=0 
      [res1,~]=size(Res_1);
      if res1~=0
         Sign_v=Res_1(:,1);
         Flag_AA=Res_1(:,2:14);
         other_info=Res_1(:,15:end);
         for i=1:res1
             F=nonzeros(Flag_AA(i,:))';
             if Sign_v(i)==510
                Hydro_AA='GAVLIMFWPYRKH';
                Flag_AA1{i}=Hydro_AA(F);
                Sign1{i}='>';
             elseif Sign_v(i)==100
                Pola_AA='GASTQNYDE';
                Flag_AA1{i}=Pola_AA(F);
                Sign1{i}='<';
             elseif Sign_v(i)==710
                Flag_AA1{i}='charge';
                Sign1{i}='>';
             end
         end
      end
      Flag=cell2table(Flag_AA1');
      Sign=cell2table(Sign1');
      Flag.Properties.VariableNames{'Var1'}='Flag';
      Sign.Properties.VariableNames{'Var1'}='Sign';
      info=array2table(other_info);
      R=cell(res1,1);
      R(:)={Region(re3)};
      Singleflag=[R,Flag,Sign,info];
   else
      Singleflag=[];
   end
   
  else
   Singleflag=[];
end
re3
  Finalflag=[Finalflag;Singleflag];
end
Region1={'H1','H1','H1','H1','H1','H1','H1','H1','H1','H1','H2','H2','H2','H2','H3'};
Region2={'H2','H2','H2','H2','H2','H2','H3','H3','H3','L1','H3','H3','H3','L1','L1'};
Region3={'H3','H3','H3','L1','L1','L2','L1','L1','L2','L2','L1','L1','L2','L2','L2'};
Region4={'L1','L2','L3','L2','L3','L3','L2','L3','L3','L3','L2','L3','L3','L3','L3'};
Region={'H1H2H3L1','H1H2H3L2','H1H2H3L3','H1H2L1L2','H1H2L1L3','H1H2L2L3','H1H3L2L3','H1H3L1L2','H1H3L2L3','H1L1L2L3','H2H3L1L2','H2H3L1L3','H2H3L2L3','H2L1L2L3','H3L1L2L3'};
for re2=1:15
    clear opt opt_0 opt_1 opt_2 Flag_AA1 Sign1 Sign_v Res_1 Row
warning('off','all');
CAB1=xlsread('Clinical_AA',Region1{re2},'A1:T200');
%HAB=xlsread('human antibody_sequence_AA comp_newlimit_v1',Region,'A2:T20013');
CP1=xlsread('Clinical_AA_PSR',Region1{re2},'A1:T300');
CA1=xlsread('Clinical_AA_ACSINS',Region1{re2},'A1:T300');
CC1=xlsread('Clinical_AA_CSI',Region1{re2},'A1:T300');
CE1=xlsread('Clinical_AA_ELISA',Region1{re2},'A1:T300');
CB1=xlsread('Clinical_AA_BVP',Region1{re2},'A1:T300');
AAB1=xlsread('Human_AA_100sim',Region1{re2},'A1:T2000');

CAB2=xlsread('Clinical_AA',Region2{re2},'A1:T200');
%HAB=xlsread('human antibody_sequence_AA comp_newlimit_v1',Region,'A2:T20013');
CP2=xlsread('Clinical_AA_PSR',Region2{re2},'A1:T300');
CA2=xlsread('Clinical_AA_ACSINS',Region2{re2},'A1:T300');
CC2=xlsread('Clinical_AA_CSI',Region2{re2},'A1:T300');
CE2=xlsread('Clinical_AA_ELISA',Region2{re2},'A1:T300');
CB2=xlsread('Clinical_AA_BVP',Region2{re2},'A1:T300');
AAB2=xlsread('Human_AA_100sim',Region2{re2},'A1:T2000');

CAB4=xlsread('Clinical_AA',Region3{re2},'A1:T200');
%HAB=xlsread('human antibody_sequence_AA comp_newlimit_v1',Region,'A2:T20013');
CP4=xlsread('Clinical_AA_PSR',Region3{re2},'A1:T300');
CA4=xlsread('Clinical_AA_ACSINS',Region3{re2},'A1:T300');
CC4=xlsread('Clinical_AA_CSI',Region3{re2},'A1:T300');
CE4=xlsread('Clinical_AA_ELISA',Region3{re2},'A1:T300');
CB4=xlsread('Clinical_AA_BVP',Region3{re2},'A1:T300');
AAB4=xlsread('Human_AA_100sim',Region3{re2},'A1:T2000');

CAB3=xlsread('Clinical_AA',Region4{re2},'A1:T200');
%HAB=xlsread('human antibody_sequence_AA comp_newlimit_v1',Region,'A2:T20013');
CP3=xlsread('Clinical_AA_PSR',Region4{re2},'A1:T300');
CA3=xlsread('Clinical_AA_ACSINS',Region4{re2},'A1:T300');
CC3=xlsread('Clinical_AA_CSI',Region4{re2},'A1:T300');
CE3=xlsread('Clinical_AA_ELISA',Region4{re2},'A1:T300');
CB3=xlsread('Clinical_AA_BVP',Region4{re2},'A1:T300');
AAB3=xlsread('Human_AA_100sim',Region4{re2},'A1:T2000');

CAB=CAB1+CAB2+CAB3+CAB4;
CP=CP1+CP2+CP3+CP4;
CA=CA1+CA2+CA3+CA4;
CC=CC1+CC2+CC3+CC4;
CE=CE1+CE2+CE3+CE4;
CB=CB1+CB2+CB3+CB4;
AAB=AAB1+AAB2+AAB3+AAB4;
adihs=[];adils=[];clis_ph=[];clis_pl=[];clis_ch=[];clis_cl=[];clis_eh=[];clis_el=[];clis_bh=[];clis_bl=[];
clis_ah=[];clis_al=[];clihs=[];clils=[];adi=[];clim=[];


adi=xlsread('2ndset5_combine_singleflags_v2_new','adimab','A2:CV100');
clis_p=xlsread('2ndset5_combine_singleflags_v2_new','PSR','A2:CV20');
clis_a=xlsread('2ndset5_combine_singleflags_v2_new','ACSINS','A2:CV20');
clis_c=xlsread('2ndset5_combine_singleflags_v2_new','CSI','A2:CV20');
clis_e=xlsread('2ndset5_combine_singleflags_v2_new','ELISA','A2:CV20');
clis_b=xlsread('2ndset5_combine_singleflags_v2_new','BVP','A2:CV20');
clim=xlsread('2ndset5_combine_singleflags_v2_new','clinical','A2:CV100');
ngp=sum(clis_p(:,1)<=gp);nbp=sum(clis_p(:,1)>gp);
nga=sum(clis_a(:,1)<=ga);nba=sum(clis_a(:,1)>ga);
ngc=sum(clis_c(:,1)<=gc);nbc=sum(clis_c(:,1)>gc);
nge=sum(clis_e(:,1)<=ge);nbe=sum(clis_e(:,1)>ge);
ngb=sum(clis_b(:,1)<=gb);nbb=sum(clis_b(:,1)>gb);
nha=sum(adi(:,1)<=ha);nla=sum(adi(:,1)>ha);
nhc=sum(clim(:,1)<=hc);nlc=sum(clim(:,1)>hc);
%hv=1033;lv=59;hr=104;lr=6;highP=0.15;
A=CAB(Total2,1);C=CAB(Total2,2);D=CAB(Total2,3);E=CAB(Total2,4);F=CAB(Total2,5);G=CAB(Total2,6);H=CAB(Total2,7);
I=CAB(Total2,8);K=CAB(Total2,9);L=CAB(Total2,10);M=CAB(Total2,11);N=CAB(Total2,12);P=CAB(Total2,13);Q=CAB(Total2,14);
R=CAB(Total2,15);S=CAB(Total2,16);T=CAB(Total2,17);V=CAB(Total2,18);W=CAB(Total2,19);Y=CAB(Total2,20);
netcharge=R+K+0.1*H-D-E;
max_nc=floor(max(netcharge));
%% Human Antibody information PSR
A1=CP(PSR2,1);C1=CP(PSR2,2);D1=CP(PSR2,3);E1=CP(PSR2,4);F1=CP(PSR2,5);G1=CP(PSR2,6);H1=CP(PSR2,7);
I1=CP(PSR2,8);K1=CP(PSR2,9);L1=CP(PSR2,10);M1=CP(PSR2,11);N1=CP(PSR2,12);P1=CP(PSR2,13);Q1=CP(PSR2,14);
R1=CP(PSR2,15);S1=CP(PSR2,16);T1=CP(PSR2,17);V1=CP(PSR2,18);W1=CP(PSR2,19);Y1=CP(PSR2,20);
Pnetcharge=R1+K1+0.1*H1-D1-E1;
%% Human Antibody information ACSINS
A2=CA(ACSINS2,1);C2=CA(ACSINS2,2);D2=CA(ACSINS2,3);E2=CA(ACSINS2,4);F2=CA(ACSINS2,5);G2=CA(ACSINS2,6);H2=CA(ACSINS2,7);
I2=CA(ACSINS2,8);K2=CA(ACSINS2,9);L2=CA(ACSINS2,10);M2=CA(ACSINS2,11);N2=CA(ACSINS2,12);P2=CA(ACSINS2,13);Q2=CA(ACSINS2,14);
R2=CA(ACSINS2,15);S2=CA(ACSINS2,16);T2=CA(ACSINS2,17);V2=CA(ACSINS2,18);W2=CA(ACSINS2,19);Y2=CA(ACSINS2,20);
Anetcharge=R2+K2+0.1*H2-D2-E2;
%% Human Antibody information CSI
A3=CC(CSI2,1);D3=CC(CSI2,3);E3=CC(CSI2,4);F3=CC(CSI2,5);G3=CC(CSI2,6);H3=CC(CSI2,7);
I3=CC(CSI2,8);K3=CC(CSI2,9);L3=CC(CSI2,10);M3=CC(CSI2,11);N3=CC(CSI2,12);P3=CC(CSI2,13);Q3=CC(CSI2,14);
R3=CC(CSI2,15);S3=CC(CSI2,16);T3=CC(CSI2,17);V3=CC(CSI2,18);W3=CC(CSI2,19);Y3=CC(CSI2,20);
Cnetcharge=R3+K3+0.1*H3-D3-E3;
%% Human Antibody information ELISA
A4=CE(ELISA2,1);C4=CE(ELISA2,2);D4=CE(ELISA2,3);E4=CE(ELISA2,4);F4=CE(ELISA2,5);G4=CE(ELISA2,6);H4=CE(ELISA2,7);
I4=CE(ELISA2,8);K4=CE(ELISA2,9);L4=CE(ELISA2,10);M4=CE(ELISA2,11);N4=CE(ELISA2,12);P4=CE(ELISA2,13);Q4=CE(ELISA2,14);
R4=CE(ELISA2,15);S4=CE(ELISA2,16);T4=CE(ELISA2,17);V4=CE(ELISA2,18);W4=CE(ELISA2,19);Y4=CE(ELISA2,20);
Enetcharge=R4+K4+0.1*H4-D4-E4;
%% Human Antibody information ELISA
A5=CB(BVP2,1);C5=CB(BVP2,2);D5=CB(BVP2,3);E5=CB(BVP2,4);F5=CB(BVP2,5);G5=CB(BVP2,6);H5=CB(BVP2,7);
I5=CB(BVP2,8);K5=CB(BVP2,9);L5=CB(BVP2,10);M5=CB(BVP2,11);N5=CB(BVP2,12);P5=CB(BVP2,13);Q5=CB(BVP2,14);
R5=CB(BVP2,15);S5=CB(BVP2,16);T5=CB(BVP2,17);V5=CB(BVP2,18);W5=CB(BVP2,19);Y5=CB(BVP2,20);
Bnetcharge=R5+K5+0.1*H5-D5-E5;
%% Human B cell from Adimab
ABA=AAB(train2,1);ABC=AAB(train2,2);ABD=AAB(train2,3);ABE=AAB(train2,4);ABF=AAB(train2,5);ABG=AAB(train2,6);ABH=AAB(train2,7);
ABI=AAB(train2,8);ABK=AAB(train2,9);ABL=AAB(train2,10);ABM=AAB(train2,11);ABN=AAB(train2,12);ABP=AAB(train2,13);ABQ=AAB(train2,14);
ABR=AAB(train2,15);ABS=AAB(train2,16);ABT=AAB(train2,17);ABV=AAB(train2,18);ABW=AAB(train2,19);ABY=AAB(train2,20);
ABnetcharge=ABR+ABK+0.1*ABH-ABD-ABE;
%% combination
% hydrophybic residue GAVLIMFWPYRKH
Hydro=[G,A,V,L,I,M,F,W,P,Y,R,K,H];
PHydro=[G1,A1,V1,L1,I1,M1,F1,W1,P1,Y1,R1,K1,H1];
AHydro=[G2,A2,V2,L2,I2,M2,F2,W2,P2,Y2,R2,K2,H2];
CHydro=[G3,A3,V3,L3,I3,M3,F3,W3,P3,Y3,R3,K3,H3];
EHydro=[G4,A4,V4,L4,I4,M4,F4,W4,P4,Y4,R4,K4,H4];
BHydro=[G5,A5,V5,L5,I5,M5,F5,W5,P5,Y5,R5,K5,H5];
AdHydro=[ABG,ABA,ABV,ABL,ABI,ABM,ABF,ABW,ABP,ABY,ABR,ABK,ABH];
Hydro_num=[1,2,3,4,5,6,7,8,9,10,11,12,13];
for i=1:13
    Hydro1{i}=combntns(1:13,i);
    [r(i),~]=size(Hydro1{i});
end    
for k=1:13
    for i=1:1:r(k)
        Hydro_1{k}(:,i)=sum(Hydro(:,Hydro1{k}(i,:)),2);
        PHydro_1{k}(:,i)=sum(PHydro(:,Hydro1{k}(i,:)),2);
        AHydro_1{k}(:,i)=sum(AHydro(:,Hydro1{k}(i,:)),2);
        CHydro_1{k}(:,i)=sum(CHydro(:,Hydro1{k}(i,:)),2);
        EHydro_1{k}(:,i)=sum(EHydro(:,Hydro1{k}(i,:)),2);
        BHydro_1{k}(:,i)=sum(BHydro(:,Hydro1{k}(i,:)),2);
        AdHydro_1{k}(:,i)=sum(AdHydro(:,Hydro1{k}(i,:)),2);
    end
    [~,r_1(k)]=size(Hydro_1{k});
    for j=1:1:r_1(k)
        maxbind_hy{k}(j)=max(Hydro_1{k}(:,j));
        minbind_hy{k}(j)=min(Hydro_1{k}(:,j));
        hy{k}(j)=(maxbind_hy{k}(j)-minbind_hy{k}(j))+1;
    end        
end    
for x=1:1:13
    q1=1;
    opt_0=[];
    for i=1:1:r(x)
        for a=minbind_hy{x}(i):1:maxbind_hy{x}(i)
            flags=Hydro_1{x}(:,i)>a;
            for f=1:50
                Fl=flags*1;
                Fl(clim(:,f),:)=[];
                Fl_1(:,f)=Fl;
            end
            averghc=(sum(flags(1:(hc),:),1)/(hc)*100);
            averglc=(sum(flags((hc+1):end,:),1)/(lc)*100);
            stdhc=std(sum(Fl_1(1:(hc-nhc)),1)/(hc-nhc)*100);
            stdlc=std(sum(Fl_1((hc-nhc+1):end),1)/(lc-nlc)*100);
            if averghc-stdhc<high && averglc+stdlc>low
               mAb=sum(flags)/(hc+lc)*100;
               highspec=averghc;
               lowspec=averglc;
               flag=a;
               Flags=flags';
               table(1,1)=sum(flags(1:hc)<1);
               table(2,1)=hc-table(1,1);
               table(1,2)=sum(flags(hc+1:end)<1);
               table(2,2)=lc-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=hc+lc;
               pvalue=TwoTailed(table);
               
               Adflags=AdHydro_1{x}(:,i)>a;
               Pflags=PHydro_1{x}(:,i)>a;
               Aflags=AHydro_1{x}(:,i)>a;
               Cflags=CHydro_1{x}(:,i)>a;
               Eflags=EHydro_1{x}(:,i)>a;
               Bflags=BHydro_1{x}(:,i)>a;
               
               table(1,1)=sum(Adflags(1:ha)<1);
               table(2,1)=ha-table(1,1);
               table(1,2)=sum(Adflags(ha+1:end)<1);
               table(2,2)=la-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=ha+la;
               pvalue1=TwoTailed(table);
               
               
               for f=1:50
                   Ad=Adflags;
                   Pflags1=Pflags;
                   Aflags1=Aflags;
                   Cflags1=Cflags;
                   Eflags1=Eflags;
                   Bflags1=Bflags;
                   
                   Ad(adi(:,f),:)=[];
                   Pflags1(clis_p(:,f),:)=[];
                   Aflags1(clis_a(:,f),:)=[];
                   Cflags1(clis_c(:,f),:)=[];
                   Eflags1(clis_e(:,f),:)=[];
                   Bflags1(clis_b(:,f),:)=[];
                   
                   Ad_1(:,f)=Ad;
                   Pflags1_1(:,f)=Pflags1;
                   Aflags1_1(:,f)=Aflags1;
                   Cflags1_1(:,f)=Cflags1;
                   Eflags1_1(:,f)=Eflags1;
                   Bflags1_1(:,f)=Bflags1;
               end
               
               AdflagH=sum(Ad_1(1:(ha-nha),:),1)/(ha-nha)*100;
               AdflagL=sum(Ad_1((ha-nha+1):end,:),1)/(la-nla)*100;
               
               Pflag1H=sum(Pflags1_1(1:(gp-ngp),:),1)/(gp-ngp)*100;
               Aflag1H=sum(Aflags1_1(1:(ga-nga),:),1)/(ga-nga)*100;
               Cflag1H=sum(Cflags1_1(1:(gc-ngc),:),1)/(gc-ngc)*100;
               Eflag1H=sum(Eflags1_1(1:(ge-nge),:),1)/(ge-nge)*100;
               Bflag1H=sum(Bflags1_1(1:(gb-ngb),:),1)/(gb-ngb)*100;
               
               Pflag1L=sum(Pflags1_1((gp-ngp+1):end,:),1)/(bp-nbp)*100;
               Aflag1L=sum(Aflags1_1((ga-nga+1):end,:),1)/(ba-nba)*100;
               Cflag1L=sum(Cflags1_1((gc-ngc+1):end,:),1)/(bc-nbc)*100;
               Eflag1L=sum(Eflags1_1((ge-nge+1):end,:),1)/(be-nbe)*100;
               Bflag1L=sum(Bflags1_1((gb-ngb+1):end,:),1)/(bb-nbb)*100;
               
               avergha=sum(Adflags(1:ha,:),1)/ha*100;
               avergla=sum(Adflags(ha+1:end,:),1)/la*100;
               averghcp=sum(Pflags(1:gp,:),1)/(gp)*100;
               averghca=sum(Aflags(1:ga,:),1)/(ga)*100;
               averghcc=sum(Cflags(1:gc,:),1)/(gc)*100;
               averghce=sum(Eflags(1:ge,:),1)/(ge)*100;
               averghcb=sum(Bflags(1:gb,:),1)/(gb)*100;
               averglcp=sum(Pflags(gp+1:end,:),1)/(bp)*100;
               averglca=sum(Aflags(ga+1:end,:),1)/(ba)*100;
               averglcc=sum(Cflags(gc+1:end,:),1)/(bc)*100;
               averglce=sum(Eflags(ge+1:end,:),1)/(be)*100;
               averglcb=sum(Bflags(gb+1:end,:),1)/(bb)*100;
              
               stdha=std(AdflagH);
               stdla=std(AdflagL);
               stdhcp=std(Pflag1H);
               stdlcp=std(Pflag1L);
               stdhca=std(Aflag1H);
               stdlca=std(Aflag1L);
               stdhcc=std(Cflag1H);
               stdlcc=std(Cflag1L);
               stdhce=std(Eflag1H);
               stdlce=std(Eflag1L);
               stdhcb=std(Bflag1H);
               stdlcb=std(Bflag1L);
               
                  
               if pvalue<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10 
                  opt{x}(q1,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
                  opt_0=zeros(1,13-x);
                  opt_1{x}(q1,:)=[510,opt_0,Hydro_num(Hydro1{x}(i,:))];
                  opt_2{x}(q1,:)=Flags;
                  q1=q1+1;
                  
               elseif pvalue1<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10  
                  opt{x}(q1,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
                  opt_0=zeros(1,13-x);
                  opt_1{x}(q1,:)=[510,opt_0,Hydro_num(Hydro1{x}(i,:))];
                  opt_2{x}(q1,:)=Flags;
                  q1=q1+1;
               end
            end
        end
    end
end
%% combination 
% polar residue GAVLIMFWP
Polar=[G,A,S,T,Q,N,Y,D,E];
PPolar=[G1,A1,S1,T1,Q1,N1,Y1,D1,E1];
APolar=[G2,A2,S2,T2,Q2,N2,Y2,D2,E2];
CPolar=[G3,A3,S3,T3,Q3,N3,Y3,D3,E3];
EPolar=[G4,A4,S4,T4,Q4,N4,Y4,D4,E4];
BPolar=[G5,A5,S5,T5,Q5,N5,Y5,D5,E5];
AdPolar=[ABG,ABA,ABS,ABT,ABQ,ABN,ABY,ABD,ABE];
Polar_num=[1,2,3,4,5,6,7,8,9];
for i=1:9
    Polar1{i}=combntns(1:9,i);
    [r_p(i),~]=size(Polar1{i});
end 
for k=1:9
    for i=1:1:r_p(k)
        Polar_1{k}(:,i)=sum(Polar(:,Polar1{k}(i,:)),2);
        PPolar_1{k}(:,i)=sum(PPolar(:,Polar1{k}(i,:)),2);
        APolar_1{k}(:,i)=sum(APolar(:,Polar1{k}(i,:)),2);
        CPolar_1{k}(:,i)=sum(CPolar(:,Polar1{k}(i,:)),2);
        EPolar_1{k}(:,i)=sum(EPolar(:,Polar1{k}(i,:)),2);
        BPolar_1{k}(:,i)=sum(BPolar(:,Polar1{k}(i,:)),2);
        AdPolar_1{k}(:,i)=sum(AdPolar(:,Polar1{k}(i,:)),2);
    end
    [~,r_p1(k)]=size(Polar_1{k});
    for j=1:1:r_p1(k)
        maxbind_po{k}(j)=max(Polar_1{k}(:,j));
        minbind_po{k}(j)=min(Polar_1{k}(:,j));
        po{k}(j)=(maxbind_po{k}(j)-minbind_po{k}(j))+1;
    end        
end

for x=1:1:9
    q2=1;
    opt_0=[];
    for i=1:1:r_p(x)
        for a=minbind_po{x}(i):1:maxbind_po{x}(i)
             flags=Polar_1{x}(:,i)<a;
            for f=1:50
                Fl=flags*1;
                Fl(clim(:,f),:)=[];
                Fl_1(:,f)=Fl;
            end
            averghc=sum(flags(1:(hc),:),1)/(hc)*100;
            averglc=(sum(flags((hc+1):end,:),1)/(lc)*100);
            stdhc=std(sum(Fl_1(1:(hc-nhc)),1)/(hc-nhc)*100);
            stdlc=std(sum(Fl_1((hc-nhc+1):end),1)/(lc-nlc)*100);
            if averghc-stdhc<high && averglc+stdlc>low
               mAb=sum(flags)/(hc+lc)*100;
               highspec=averghc;
               lowspec=averglc;
               flag=a;
               Flags=flags';
               table(1,1)=sum(flags(1:hc)<1);
               table(2,1)=hc-table(1,1);
               table(1,2)=sum(flags(hc+1:end)<1);
               table(2,2)=lc-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=hc+lc;
               pvalue=TwoTailed(table);
               
               Adflags=AdPolar_1{x}(:,i)<a;
               Pflags=PPolar_1{x}(:,i)<a;
               Aflags=APolar_1{x}(:,i)<a;
               Cflags=CPolar_1{x}(:,i)<a;
               Eflags=EPolar_1{x}(:,i)<a;
               Bflags=BPolar_1{x}(:,i)<a;
               
               table(1,1)=sum(Adflags(1:ha)<1);
               table(2,1)=ha-table(1,1);
               table(1,2)=sum(Adflags(ha+1:end)<1);
               table(2,2)=la-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=ha+la;
               pvalue1=TwoTailed(table);
               
               for f=1:50
                   Ad=Adflags;
                   Pflags1=Pflags;
                   Aflags1=Aflags;
                   Cflags1=Cflags;
                   Eflags1=Eflags;
                   Bflags1=Bflags;
                   
                   Ad(adi(:,f),:)=[];
                   Pflags1(clis_p(:,f),:)=[];
                   Aflags1(clis_a(:,f),:)=[];
                   Cflags1(clis_c(:,f),:)=[];
                   Eflags1(clis_e(:,f),:)=[];
                   Bflags1(clis_b(:,f),:)=[];
                   
                   Ad_1(:,f)=Ad;
                   Pflags1_1(:,f)=Pflags1;
                   Aflags1_1(:,f)=Aflags1;
                   Cflags1_1(:,f)=Cflags1;
                   Eflags1_1(:,f)=Eflags1;
                   Bflags1_1(:,f)=Bflags1;
               end
               
               AdflagH=sum(Ad_1(1:(ha-nha),:),1)/(ha-nha)*100;
               AdflagL=sum(Ad_1((ha-nha+1):end,:),1)/(la-nla)*100;
               
               Pflag1H=sum(Pflags1_1(1:(gp-ngp),:),1)/(gp-ngp)*100;
               Aflag1H=sum(Aflags1_1(1:(ga-nga),:),1)/(ga-nga)*100;
               Cflag1H=sum(Cflags1_1(1:(gc-ngc),:),1)/(gc-ngc)*100;
               Eflag1H=sum(Eflags1_1(1:(ge-nge),:),1)/(ge-nge)*100;
               Bflag1H=sum(Bflags1_1(1:(gb-ngb),:),1)/(gb-ngb)*100;
               
               Pflag1L=sum(Pflags1_1((gp-ngp+1):end,:),1)/(bp-nbp)*100;
               Aflag1L=sum(Aflags1_1((ga-nga+1):end,:),1)/(ba-nba)*100;
               Cflag1L=sum(Cflags1_1((gc-ngc+1):end,:),1)/(bc-nbc)*100;
               Eflag1L=sum(Eflags1_1((ge-nge+1):end,:),1)/(be-nbe)*100;
               Bflag1L=sum(Bflags1_1((gb-ngb+1):end,:),1)/(bb-nbb)*100;
               
               avergha=sum(Adflags(1:ha,:),1)/ha*100;
               avergla=sum(Adflags(ha+1:end,:),1)/la*100;
               averghcp=sum(Pflags(1:gp,:),1)/(gp)*100;
               averghca=sum(Aflags(1:ga,:),1)/(ga)*100;
               averghcc=sum(Cflags(1:gc,:),1)/(gc)*100;
               averghce=sum(Eflags(1:ge,:),1)/(ge)*100;
               averghcb=sum(Bflags(1:gb,:),1)/(gb)*100;
               averglcp=sum(Pflags(gp+1:end,:),1)/(bp)*100;
               averglca=sum(Aflags(ga+1:end,:),1)/(ba)*100;
               averglcc=sum(Cflags(gc+1:end,:),1)/(bc)*100;
               averglce=sum(Eflags(ge+1:end,:),1)/(be)*100;
               averglcb=sum(Bflags(gb+1:end,:),1)/(bb)*100;
              
               stdha=std(AdflagH);
               stdla=std(AdflagL);
               stdhcp=std(Pflag1H);
               stdlcp=std(Pflag1L);
               stdhca=std(Aflag1H);
               stdlca=std(Aflag1L);
               stdhcc=std(Cflag1H);
               stdlcc=std(Cflag1L);
               stdhce=std(Eflag1H);
               stdlce=std(Eflag1L);
               stdhcb=std(Bflag1H);
               stdlcb=std(Bflag1L);  
               if pvalue<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10 
                  opt{x+13+13}(q2,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
                  opt_0=zeros(1,13-x);
                  opt_1{x+13+13}(q2,:)=[100,opt_0,Polar_num(Polar1{x}(i,:))];
                  opt_2{x+13+13}(q2,:)=Flags;
                  q2=q2+1;
               elseif pvalue1<0.05 &&averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10  
                  opt{x+13+13}(q2,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
                  opt_0=zeros(1,13-x);
                  opt_1{x+13+13}(q2,:)=[100,opt_0,Polar_num(Polar1{x}(i,:))];
                  opt_2{x+13+13}(q2,:)=Flags;
                  q2=q2+1;
               end
            end
        end
    end
end

%% single amino acid
Single=netcharge;
PSingle=Pnetcharge;
ASingle=Anetcharge;
CSingle=Cnetcharge;
ESingle=Enetcharge;
BSingle=Bnetcharge;
AdSingle=ABnetcharge;
y=1;
for a=0:1:max_nc
    flags=Single>a;
    for f=1:50
        Fl=flags*1;
        Fl(clim(:,f),:)=[];
        Fl_1(:,f)=Fl;
    end
    averghc=(sum(flags(1:(hc),:),1)/(hc)*100);
            averglc=(sum(flags((hc+1):end,:),1)/(lc)*100);
            stdhc=std(sum(Fl_1(1:(hc-nhc)),1)/(hc-nhc)*100);
            stdlc=std(sum(Fl_1((hc-nhc+1):end),1)/(lc-nlc)*100);
            if averghc-stdhc<high && averglc+stdlc>low
               mAb=sum(flags)/(hc+lc)*100;
               highspec=averghc;
               lowspec=averglc;
               flag=a;
               Flags=flags';
               table(1,1)=sum(flags(1:hc)<1);
               table(2,1)=hc-table(1,1);
               table(1,2)=sum(flags(hc+1:end)<1);
               table(2,2)=lc-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=hc+lc;
               pvalue=TwoTailed(table);
 
       Adflags=AdSingle>a;
       Pflags=PSingle>a;
       Aflags=ASingle>a;
       Cflags=CSingle>a;
       Eflags=ESingle>a;
       Bflags=BSingle>a;
       
       table(1,1)=sum(Adflags(1:ha)<1);
               table(2,1)=ha-table(1,1);
               table(1,2)=sum(Adflags(ha+1:end)<1);
               table(2,2)=la-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=ha+la;
               pvalue1=TwoTailed(table);
       for f=1:50
                   Ad=Adflags;
                   Pflags1=Pflags;
                   Aflags1=Aflags;
                   Cflags1=Cflags;
                   Eflags1=Eflags;
                   Bflags1=Bflags;
                   
                   Ad(adi(:,f),:)=[];
                   Pflags1(clis_p(:,f),:)=[];
                   Aflags1(clis_a(:,f),:)=[];
                   Cflags1(clis_c(:,f),:)=[];
                   Eflags1(clis_e(:,f),:)=[];
                   Bflags1(clis_b(:,f),:)=[];
                   
                   Ad_1(:,f)=Ad;
                   Pflags1_1(:,f)=Pflags1;
                   Aflags1_1(:,f)=Aflags1;
                   Cflags1_1(:,f)=Cflags1;
                   Eflags1_1(:,f)=Eflags1;
                   Bflags1_1(:,f)=Bflags1;
               end
               
               AdflagH=sum(Ad_1(1:(ha-nha),:),1)/(ha-nha)*100;
               AdflagL=sum(Ad_1((ha-nha+1):end,:),1)/(la-nla)*100;
               
               Pflag1H=sum(Pflags1_1(1:(gp-ngp),:),1)/(gp-ngp)*100;
               Aflag1H=sum(Aflags1_1(1:(ga-nga),:),1)/(ga-nga)*100;
               Cflag1H=sum(Cflags1_1(1:(gc-ngc),:),1)/(gc-ngc)*100;
               Eflag1H=sum(Eflags1_1(1:(ge-nge),:),1)/(ge-nge)*100;
               Bflag1H=sum(Bflags1_1(1:(gb-ngb),:),1)/(gb-ngb)*100;
               
               Pflag1L=sum(Pflags1_1((gp-ngp+1):end,:),1)/(bp-nbp)*100;
               Aflag1L=sum(Aflags1_1((ga-nga+1):end,:),1)/(ba-nba)*100;
               Cflag1L=sum(Cflags1_1((gc-ngc+1):end,:),1)/(bc-nbc)*100;
               Eflag1L=sum(Eflags1_1((ge-nge+1):end,:),1)/(be-nbe)*100;
               Bflag1L=sum(Bflags1_1((gb-ngb+1):end,:),1)/(bb-nbb)*100;
               
               avergha=sum(Adflags(1:ha,:),1)/ha*100;
               avergla=sum(Adflags(ha+1:end,:),1)/la*100;
               averghcp=sum(Pflags(1:gp,:),1)/(gp)*100;
               averghca=sum(Aflags(1:ga,:),1)/(ga)*100;
               averghcc=sum(Cflags(1:gc,:),1)/(gc)*100;
               averghce=sum(Eflags(1:ge,:),1)/(ge)*100;
               averghcb=sum(Bflags(1:gb,:),1)/(gb)*100;
               averglcp=sum(Pflags(gp+1:end,:),1)/(bp)*100;
               averglca=sum(Aflags(ga+1:end,:),1)/(ba)*100;
               averglcc=sum(Cflags(gc+1:end,:),1)/(bc)*100;
               averglce=sum(Eflags(ge+1:end,:),1)/(be)*100;
               averglcb=sum(Bflags(gb+1:end,:),1)/(bb)*100;
              
               stdha=std(AdflagH);
               stdla=std(AdflagL);
               stdhcp=std(Pflag1H);
               stdlcp=std(Pflag1L);
               stdhca=std(Aflag1H);
               stdlca=std(Aflag1L);
               stdhcc=std(Cflag1H);
               stdlcc=std(Cflag1L);
               stdhce=std(Eflag1H);
               stdlce=std(Eflag1L);
               stdhcb=std(Bflag1H);
               stdlcb=std(Bflag1L);
       if pvalue<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10  
          opt{45}(y,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
          opt_1{45}(y,:)=[710,0,0,0,0,0,0,0,0,0,0,0,0,1];
          opt_2{45}(y,:)=[Flags];
          y=y+1;
       elseif pvalue1<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10  
          opt{45}(y,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
          opt_1{45}(y,:)=[710,0,0,0,0,0,0,0,0,0,0,0,0,1];
          opt_2{45}(y,:)=[Flags];
          y=y+1;
       end
    end
end
               
exist opt               
if ans~=0
   [~,r]=size(opt);
   Res=[];
   for i=1:r
       Res=[Res;opt_1{i} opt{i} opt_2{i}];
   end
   Res(isnan(Res(:,1)),:) = [] ;

   [e,~]=size(Res);
   Fx=Res(:,1:14);
   Px=Res(:,15:19);
   q=0;
   for i=1:e-1
       for j=i+1:e
           p=ismember(Fx(i,:),Fx(j,:));
           p1=sum(p);
           if p1==14 && Px(i,1)==Px(j,1)&& Px(i,2)==Px(j,2) && Px(i,3)==Px(j,3)
              q=q+1;
              Row(q)=j;
           end
       end
   end

   if q~=0
      Row=unique(Row);
      Res(Row,:) = []; 
      Res(:,18)=-1*Res(:,18);
      Res(:,20)=-1*Res(:,20);
      Res(:,22)=-1*Res(:,22);
      Res(:,24)=-1*Res(:,24);
      Res(:,26)=-1*Res(:,26);
      Res(:,28)=-1*Res(:,28);
      Res(:,30)=-1*Res(:,30);
      Res=sortrows(Res,[16,18,30,20,22,24,26,28]);
      Res(:,18)=-1*Res(:,18);
      Res(:,20)=-1*Res(:,20);
      Res(:,22)=-1*Res(:,22);
      Res(:,24)=-1*Res(:,24);
      Res(:,26)=-1*Res(:,26);
      Res(:,28)=-1*Res(:,28);
      Res(:,30)=-1*Res(:,30);
      [res,~]=size(Res);
      Sig=zeros(res,1);
      for i=1:res
          pc=i/res*0.05;
          if Res(i,16)<pc
             Sig(i,:)=1;
          end
      end
      s=find(Sig==1);
      [rs,~]=size(s);
      if rs==0
         disp('No results')
      else
      max1=max(s);
      Res_1=Res(1:max1,:); % results from FDR
      end
   else
      Res(:,18)=-1*Res(:,18);
      Res(:,20)=-1*Res(:,20);
      Res(:,22)=-1*Res(:,22);
      Res(:,24)=-1*Res(:,24);
      Res(:,26)=-1*Res(:,26);
      Res(:,28)=-1*Res(:,28);
      Res(:,30)=-1*Res(:,30);
      Res=sortrows(Res,[16,18,30,20,22,24,26,28]);
      Res(:,18)=-1*Res(:,18);
      Res(:,20)=-1*Res(:,20);
      Res(:,22)=-1*Res(:,22);
      Res(:,24)=-1*Res(:,24);
      Res(:,26)=-1*Res(:,26);
      Res(:,28)=-1*Res(:,28);
      Res(:,30)=-1*Res(:,30);
      [res,~]=size(Res);
      Sig=zeros(res,1);
      for i=1:res
          pc=i/res*0.05;
          if Res(i,16)<pc
          Sig(i,:)=1;
          end
      end
      s=find(Sig==1);
      [rs,~]=size(s);
      if rs==0
         disp('No results')
      else
         max1=max(s);
         Res_1=Res(1:max1,:); % results from FDR
      end
   end
   exist Res_1
   if ans ~=0 
      [res1,~]=size(Res_1);
      if res1~=0
         Sign_v=Res_1(:,1);
         Flag_AA=Res_1(:,2:14);
         other_info=Res_1(:,15:end);
         for i=1:res1
             F=nonzeros(Flag_AA(i,:))';
             if Sign_v(i)==510
                Hydro_AA='GAVLIMFWPYRKH';
                Flag_AA1{i}=Hydro_AA(F);
                Sign1{i}='>';
             elseif Sign_v(i)==100
                Pola_AA='GASTQNYDE';
                Flag_AA1{i}=Pola_AA(F);
                Sign1{i}='<';
             elseif Sign_v(i)==710
                Flag_AA1{i}='charge';
                Sign1{i}='>';
             end
         end
      end
      Flag=cell2table(Flag_AA1');
      Sign=cell2table(Sign1');
      Flag.Properties.VariableNames{'Var1'}='Flag';
      Sign.Properties.VariableNames{'Var1'}='Sign';
      info=array2table(other_info);
      R=cell(res1,1);
      R(:)={Region(re2)};
      Singleflag=[R,Flag,Sign,info];
   else
      Singleflag=[];
   end
   
  else
   Singleflag=[];
end
re2
  Finalflag=[Finalflag;Singleflag];
end
Region1={'H1','H1','H1','H1','H1','H2'};
Region2={'H2','H2','H2','H2','H3','H3'};
Region3={'H3','H3','H3','L1','L1','L1'};
Region4={'L1','L1','L3','L2','L2','L2'};
Region5={'L2','L3','L2','L3','L3','L3'};
Region={'H1H2H3L1L2','H1H2H3L1L3','H1H2H3L3L2','H1H2L1L2L3','H1H3L1L2L3','H2H3L1L2L3'};
for re2=1:6
    clear opt opt_0 opt_1 opt_2 Flag_AA1 Sign1 Sign_v Res_1 Row
warning('off','all');
CAB1=xlsread('Clinical_AA',Region1{re2},'A1:T200');
%HAB=xlsread('human antibody_sequence_AA comp_newlimit_v1',Region,'A2:T20013');
CP1=xlsread('Clinical_AA_PSR',Region1{re2},'A1:T300');
CA1=xlsread('Clinical_AA_ACSINS',Region1{re2},'A1:T300');
CC1=xlsread('Clinical_AA_CSI',Region1{re2},'A1:T300');
CE1=xlsread('Clinical_AA_ELISA',Region1{re2},'A1:T300');
CB1=xlsread('Clinical_AA_BVP',Region1{re2},'A1:T300');
AAB1=xlsread('Human_AA_100sim',Region1{re2},'A1:T2000');

CAB2=xlsread('Clinical_AA',Region2{re2},'A1:T200');
%HAB=xlsread('human antibody_sequence_AA comp_newlimit_v1',Region,'A2:T20013');
CP2=xlsread('Clinical_AA_PSR',Region2{re2},'A1:T300');
CA2=xlsread('Clinical_AA_ACSINS',Region2{re2},'A1:T300');
CC2=xlsread('Clinical_AA_CSI',Region2{re2},'A1:T300');
CE2=xlsread('Clinical_AA_ELISA',Region2{re2},'A1:T300');
CB2=xlsread('Clinical_AA_BVP',Region2{re2},'A1:T300');
AAB2=xlsread('Human_AA_100sim',Region2{re2},'A1:T2000');

CAB4=xlsread('Clinical_AA',Region3{re2},'A1:T200');
%HAB=xlsread('human antibody_sequence_AA comp_newlimit_v1',Region,'A2:T20013');
CP4=xlsread('Clinical_AA_PSR',Region3{re2},'A1:T300');
CA4=xlsread('Clinical_AA_ACSINS',Region3{re2},'A1:T300');
CC4=xlsread('Clinical_AA_CSI',Region3{re2},'A1:T300');
CE4=xlsread('Clinical_AA_ELISA',Region3{re2},'A1:T300');
CB4=xlsread('Clinical_AA_BVP',Region3{re2},'A1:T300');
AAB4=xlsread('Human_AA_100sim',Region3{re2},'A1:T2000');

CAB3=xlsread('Clinical_AA',Region4{re2},'A1:T200');
%HAB=xlsread('human antibody_sequence_AA comp_newlimit_v1',Region,'A2:T20013');
CP3=xlsread('Clinical_AA_PSR',Region4{re2},'A1:T300');
CA3=xlsread('Clinical_AA_ACSINS',Region4{re2},'A1:T300');
CC3=xlsread('Clinical_AA_CSI',Region4{re2},'A1:T300');
CE3=xlsread('Clinical_AA_ELISA',Region4{re2},'A1:T300');
CB3=xlsread('Clinical_AA_BVP',Region4{re2},'A1:T300');
AAB3=xlsread('Human_AA_100sim',Region4{re2},'A1:T2000');

CAB5=xlsread('Clinical_AA',Region5{re2},'A1:T200');
%HAB=xlsread('human antibody_sequence_AA comp_newlimit_v1',Region,'A2:T20013');
CP5=xlsread('Clinical_AA_PSR',Region5{re2},'A1:T300');
CA5=xlsread('Clinical_AA_ACSINS',Region5{re2},'A1:T300');
CC5=xlsread('Clinical_AA_CSI',Region5{re2},'A1:T300');
CE5=xlsread('Clinical_AA_ELISA',Region5{re2},'A1:T300');
CB5=xlsread('Clinical_AA_BVP',Region5{re2},'A1:T300');
AAB5=xlsread('Human_AA_100sim',Region5{re2},'A1:T2000');

CAB=CAB1+CAB2+CAB3+CAB4+CB5;
CP=CP1+CP2+CP3+CP4+CP5;
CA=CA1+CA2+CA3+CA4+CA5;
CC=CC1+CC2+CC3+CC4+CC5;
CE=CE1+CE2+CE3+CE4+CE5;
CB=CB1+CB2+CB3+CB4+CB5;
AAB=AAB1+AAB2+AAB3+AAB4+AAB5;
adihs=[];adils=[];clis_ph=[];clis_pl=[];clis_ch=[];clis_cl=[];clis_eh=[];clis_el=[];clis_bh=[];clis_bl=[];
clis_ah=[];clis_al=[];clihs=[];clils=[];adi=[];clim=[];


adi=xlsread('2ndset5_combine_singleflags_v2_new','adimab','A2:CV100');
clis_p=xlsread('2ndset5_combine_singleflags_v2_new','PSR','A2:CV20');
clis_a=xlsread('2ndset5_combine_singleflags_v2_new','ACSINS','A2:CV20');
clis_c=xlsread('2ndset5_combine_singleflags_v2_new','CSI','A2:CV20');
clis_e=xlsread('2ndset5_combine_singleflags_v2_new','ELISA','A2:CV20');
clis_b=xlsread('2ndset5_combine_singleflags_v2_new','BVP','A2:CV20');
clim=xlsread('2ndset5_combine_singleflags_v2_new','clinical','A2:CV100');
ngp=sum(clis_p(:,1)<=gp);nbp=sum(clis_p(:,1)>gp);
nga=sum(clis_a(:,1)<=ga);nba=sum(clis_a(:,1)>ga);
ngc=sum(clis_c(:,1)<=gc);nbc=sum(clis_c(:,1)>gc);
nge=sum(clis_e(:,1)<=ge);nbe=sum(clis_e(:,1)>ge);
ngb=sum(clis_b(:,1)<=gb);nbb=sum(clis_b(:,1)>gb);
nha=sum(adi(:,1)<=ha);nla=sum(adi(:,1)>ha);
nhc=sum(clim(:,1)<=hc);nlc=sum(clim(:,1)>hc);
%hv=1033;lv=59;hr=104;lr=6;highP=0.15;
A=CAB(Total2,1);C=CAB(Total2,2);D=CAB(Total2,3);E=CAB(Total2,4);F=CAB(Total2,5);G=CAB(Total2,6);H=CAB(Total2,7);
I=CAB(Total2,8);K=CAB(Total2,9);L=CAB(Total2,10);M=CAB(Total2,11);N=CAB(Total2,12);P=CAB(Total2,13);Q=CAB(Total2,14);
R=CAB(Total2,15);S=CAB(Total2,16);T=CAB(Total2,17);V=CAB(Total2,18);W=CAB(Total2,19);Y=CAB(Total2,20);
netcharge=R+K+0.1*H-D-E;
max_nc=floor(max(netcharge));
%% Human Antibody information PSR
A1=CP(PSR2,1);C1=CP(PSR2,2);D1=CP(PSR2,3);E1=CP(PSR2,4);F1=CP(PSR2,5);G1=CP(PSR2,6);H1=CP(PSR2,7);
I1=CP(PSR2,8);K1=CP(PSR2,9);L1=CP(PSR2,10);M1=CP(PSR2,11);N1=CP(PSR2,12);P1=CP(PSR2,13);Q1=CP(PSR2,14);
R1=CP(PSR2,15);S1=CP(PSR2,16);T1=CP(PSR2,17);V1=CP(PSR2,18);W1=CP(PSR2,19);Y1=CP(PSR2,20);
Pnetcharge=R1+K1+0.1*H1-D1-E1;
%% Human Antibody information ACSINS
A2=CA(ACSINS2,1);C2=CA(ACSINS2,2);D2=CA(ACSINS2,3);E2=CA(ACSINS2,4);F2=CA(ACSINS2,5);G2=CA(ACSINS2,6);H2=CA(ACSINS2,7);
I2=CA(ACSINS2,8);K2=CA(ACSINS2,9);L2=CA(ACSINS2,10);M2=CA(ACSINS2,11);N2=CA(ACSINS2,12);P2=CA(ACSINS2,13);Q2=CA(ACSINS2,14);
R2=CA(ACSINS2,15);S2=CA(ACSINS2,16);T2=CA(ACSINS2,17);V2=CA(ACSINS2,18);W2=CA(ACSINS2,19);Y2=CA(ACSINS2,20);
Anetcharge=R2+K2+0.1*H2-D2-E2;
%% Human Antibody information CSI
A3=CC(CSI2,1);D3=CC(CSI2,3);E3=CC(CSI2,4);F3=CC(CSI2,5);G3=CC(CSI2,6);H3=CC(CSI2,7);
I3=CC(CSI2,8);K3=CC(CSI2,9);L3=CC(CSI2,10);M3=CC(CSI2,11);N3=CC(CSI2,12);P3=CC(CSI2,13);Q3=CC(CSI2,14);
R3=CC(CSI2,15);S3=CC(CSI2,16);T3=CC(CSI2,17);V3=CC(CSI2,18);W3=CC(CSI2,19);Y3=CC(CSI2,20);
Cnetcharge=R3+K3+0.1*H3-D3-E3;
%% Human Antibody information ELISA
A4=CE(ELISA2,1);C4=CE(ELISA2,2);D4=CE(ELISA2,3);E4=CE(ELISA2,4);F4=CE(ELISA2,5);G4=CE(ELISA2,6);H4=CE(ELISA2,7);
I4=CE(ELISA2,8);K4=CE(ELISA2,9);L4=CE(ELISA2,10);M4=CE(ELISA2,11);N4=CE(ELISA2,12);P4=CE(ELISA2,13);Q4=CE(ELISA2,14);
R4=CE(ELISA2,15);S4=CE(ELISA2,16);T4=CE(ELISA2,17);V4=CE(ELISA2,18);W4=CE(ELISA2,19);Y4=CE(ELISA2,20);
Enetcharge=R4+K4+0.1*H4-D4-E4;
%% Human Antibody information ELISA
A5=CB(BVP2,1);C5=CB(BVP2,2);D5=CB(BVP2,3);E5=CB(BVP2,4);F5=CB(BVP2,5);G5=CB(BVP2,6);H5=CB(BVP2,7);
I5=CB(BVP2,8);K5=CB(BVP2,9);L5=CB(BVP2,10);M5=CB(BVP2,11);N5=CB(BVP2,12);P5=CB(BVP2,13);Q5=CB(BVP2,14);
R5=CB(BVP2,15);S5=CB(BVP2,16);T5=CB(BVP2,17);V5=CB(BVP2,18);W5=CB(BVP2,19);Y5=CB(BVP2,20);
Bnetcharge=R5+K5+0.1*H5-D5-E5;
%% Human B cell from Adimab
ABA=AAB(train2,1);ABC=AAB(train2,2);ABD=AAB(train2,3);ABE=AAB(train2,4);ABF=AAB(train2,5);ABG=AAB(train2,6);ABH=AAB(train2,7);
ABI=AAB(train2,8);ABK=AAB(train2,9);ABL=AAB(train2,10);ABM=AAB(train2,11);ABN=AAB(train2,12);ABP=AAB(train2,13);ABQ=AAB(train2,14);
ABR=AAB(train2,15);ABS=AAB(train2,16);ABT=AAB(train2,17);ABV=AAB(train2,18);ABW=AAB(train2,19);ABY=AAB(train2,20);
ABnetcharge=ABR+ABK+0.1*ABH-ABD-ABE;
%% combination
% hydrophybic residue GAVLIMFWPYRKH
Hydro=[G,A,V,L,I,M,F,W,P,Y,R,K,H];
PHydro=[G1,A1,V1,L1,I1,M1,F1,W1,P1,Y1,R1,K1,H1];
AHydro=[G2,A2,V2,L2,I2,M2,F2,W2,P2,Y2,R2,K2,H2];
CHydro=[G3,A3,V3,L3,I3,M3,F3,W3,P3,Y3,R3,K3,H3];
EHydro=[G4,A4,V4,L4,I4,M4,F4,W4,P4,Y4,R4,K4,H4];
BHydro=[G5,A5,V5,L5,I5,M5,F5,W5,P5,Y5,R5,K5,H5];
AdHydro=[ABG,ABA,ABV,ABL,ABI,ABM,ABF,ABW,ABP,ABY,ABR,ABK,ABH];
Hydro_num=[1,2,3,4,5,6,7,8,9,10,11,12,13];
for i=1:13
    Hydro1{i}=combntns(1:13,i);
    [r(i),~]=size(Hydro1{i});
end    
for k=1:13
    for i=1:1:r(k)
        Hydro_1{k}(:,i)=sum(Hydro(:,Hydro1{k}(i,:)),2);
        PHydro_1{k}(:,i)=sum(PHydro(:,Hydro1{k}(i,:)),2);
        AHydro_1{k}(:,i)=sum(AHydro(:,Hydro1{k}(i,:)),2);
        CHydro_1{k}(:,i)=sum(CHydro(:,Hydro1{k}(i,:)),2);
        EHydro_1{k}(:,i)=sum(EHydro(:,Hydro1{k}(i,:)),2);
        BHydro_1{k}(:,i)=sum(BHydro(:,Hydro1{k}(i,:)),2);
        AdHydro_1{k}(:,i)=sum(AdHydro(:,Hydro1{k}(i,:)),2);
    end
    [~,r_1(k)]=size(Hydro_1{k});
    for j=1:1:r_1(k)
        maxbind_hy{k}(j)=max(Hydro_1{k}(:,j));
        minbind_hy{k}(j)=min(Hydro_1{k}(:,j));
        hy{k}(j)=(maxbind_hy{k}(j)-minbind_hy{k}(j))+1;
    end        
end    
for x=1:1:13
    q1=1;
    opt_0=[];
    for i=1:1:r(x)
        for a=minbind_hy{x}(i):1:maxbind_hy{x}(i)
            flags=Hydro_1{x}(:,i)>a;
            for f=1:50
                Fl=flags*1;
                Fl(clim(:,f),:)=[];
                Fl_1(:,f)=Fl;
            end
            averghc=(sum(flags(1:(hc),:),1)/(hc)*100);
            averglc=(sum(flags((hc+1):end,:),1)/(lc)*100);
            stdhc=std(sum(Fl_1(1:(hc-nhc)),1)/(hc-nhc)*100);
            stdlc=std(sum(Fl_1((hc-nhc+1):end),1)/(lc-nlc)*100);
            if averghc-stdhc<high && averglc+stdlc>low
               mAb=sum(flags)/(hc+lc)*100;
               highspec=averghc;
               lowspec=averglc;
               flag=a;
               Flags=flags';
               table(1,1)=sum(flags(1:hc)<1);
               table(2,1)=hc-table(1,1);
               table(1,2)=sum(flags(hc+1:end)<1);
               table(2,2)=lc-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=hc+lc;
               pvalue=TwoTailed(table);
               
               Adflags=AdHydro_1{x}(:,i)>a;
               Pflags=PHydro_1{x}(:,i)>a;
               Aflags=AHydro_1{x}(:,i)>a;
               Cflags=CHydro_1{x}(:,i)>a;
               Eflags=EHydro_1{x}(:,i)>a;
               Bflags=BHydro_1{x}(:,i)>a;
               
               table(1,1)=sum(Adflags(1:ha)<1);
               table(2,1)=ha-table(1,1);
               table(1,2)=sum(Adflags(ha+1:end)<1);
               table(2,2)=la-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=ha+la;
               pvalue1=TwoTailed(table);
               
               
               for f=1:50
                   Ad=Adflags;
                   Pflags1=Pflags;
                   Aflags1=Aflags;
                   Cflags1=Cflags;
                   Eflags1=Eflags;
                   Bflags1=Bflags;
                   
                   Ad(adi(:,f),:)=[];
                   Pflags1(clis_p(:,f),:)=[];
                   Aflags1(clis_a(:,f),:)=[];
                   Cflags1(clis_c(:,f),:)=[];
                   Eflags1(clis_e(:,f),:)=[];
                   Bflags1(clis_b(:,f),:)=[];
                   
                   Ad_1(:,f)=Ad;
                   Pflags1_1(:,f)=Pflags1;
                   Aflags1_1(:,f)=Aflags1;
                   Cflags1_1(:,f)=Cflags1;
                   Eflags1_1(:,f)=Eflags1;
                   Bflags1_1(:,f)=Bflags1;
               end
               
               AdflagH=sum(Ad_1(1:(ha-nha),:),1)/(ha-nha)*100;
               AdflagL=sum(Ad_1((ha-nha+1):end,:),1)/(la-nla)*100;
               
               Pflag1H=sum(Pflags1_1(1:(gp-ngp),:),1)/(gp-ngp)*100;
               Aflag1H=sum(Aflags1_1(1:(ga-nga),:),1)/(ga-nga)*100;
               Cflag1H=sum(Cflags1_1(1:(gc-ngc),:),1)/(gc-ngc)*100;
               Eflag1H=sum(Eflags1_1(1:(ge-nge),:),1)/(ge-nge)*100;
               Bflag1H=sum(Bflags1_1(1:(gb-ngb),:),1)/(gb-ngb)*100;
               
               Pflag1L=sum(Pflags1_1((gp-ngp+1):end,:),1)/(bp-nbp)*100;
               Aflag1L=sum(Aflags1_1((ga-nga+1):end,:),1)/(ba-nba)*100;
               Cflag1L=sum(Cflags1_1((gc-ngc+1):end,:),1)/(bc-nbc)*100;
               Eflag1L=sum(Eflags1_1((ge-nge+1):end,:),1)/(be-nbe)*100;
               Bflag1L=sum(Bflags1_1((gb-ngb+1):end,:),1)/(bb-nbb)*100;
               
              avergha=sum(Adflags(1:ha,:),1)/ha*100;
               avergla=sum(Adflags(ha+1:end,:),1)/la*100;
               averghcp=sum(Pflags(1:gp,:),1)/(gp)*100;
               averghca=sum(Aflags(1:ga,:),1)/(ga)*100;
               averghcc=sum(Cflags(1:gc,:),1)/(gc)*100;
               averghce=sum(Eflags(1:ge,:),1)/(ge)*100;
               averghcb=sum(Bflags(1:gb,:),1)/(gb)*100;
               averglcp=sum(Pflags(gp+1:end,:),1)/(bp)*100;
               averglca=sum(Aflags(ga+1:end,:),1)/(ba)*100;
               averglcc=sum(Cflags(gc+1:end,:),1)/(bc)*100;
               averglce=sum(Eflags(ge+1:end,:),1)/(be)*100;
               averglcb=sum(Bflags(gb+1:end,:),1)/(bb)*100;
              
               stdha=std(AdflagH);
               stdla=std(AdflagL);
               stdhcp=std(Pflag1H);
               stdlcp=std(Pflag1L);
               stdhca=std(Aflag1H);
               stdlca=std(Aflag1L);
               stdhcc=std(Cflag1H);
               stdlcc=std(Cflag1L);
               stdhce=std(Eflag1H);
               stdlce=std(Eflag1L);
               stdhcb=std(Bflag1H);
               stdlcb=std(Bflag1L);
               
                  
               if pvalue<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10 
                  opt{x}(q1,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
                  opt_0=zeros(1,13-x);
                  opt_1{x}(q1,:)=[510,opt_0,Hydro_num(Hydro1{x}(i,:))];
                  opt_2{x}(q1,:)=Flags;
                  q1=q1+1;
                  
               elseif pvalue1<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10  
                  opt{x}(q1,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
                  opt_0=zeros(1,13-x);
                  opt_1{x}(q1,:)=[510,opt_0,Hydro_num(Hydro1{x}(i,:))];
                  opt_2{x}(q1,:)=Flags;
                  q1=q1+1;
               end
            end
        end
    end
end

%% combination 
% polar residue GAVLIMFWP
Polar=[G,A,S,T,Q,N,Y,D,E];
PPolar=[G1,A1,S1,T1,Q1,N1,Y1,D1,E1];
APolar=[G2,A2,S2,T2,Q2,N2,Y2,D2,E2];
CPolar=[G3,A3,S3,T3,Q3,N3,Y3,D3,E3];
EPolar=[G4,A4,S4,T4,Q4,N4,Y4,D4,E4];
BPolar=[G5,A5,S5,T5,Q5,N5,Y5,D5,E5];
AdPolar=[ABG,ABA,ABS,ABT,ABQ,ABN,ABY,ABD,ABE];
Polar_num=[1,2,3,4,5,6,7,8,9];
for i=1:9
    Polar1{i}=combntns(1:9,i);
    [r_p(i),~]=size(Polar1{i});
end 
for k=1:9
    for i=1:1:r_p(k)
        Polar_1{k}(:,i)=sum(Polar(:,Polar1{k}(i,:)),2);
        PPolar_1{k}(:,i)=sum(PPolar(:,Polar1{k}(i,:)),2);
        APolar_1{k}(:,i)=sum(APolar(:,Polar1{k}(i,:)),2);
        CPolar_1{k}(:,i)=sum(CPolar(:,Polar1{k}(i,:)),2);
        EPolar_1{k}(:,i)=sum(EPolar(:,Polar1{k}(i,:)),2);
        BPolar_1{k}(:,i)=sum(BPolar(:,Polar1{k}(i,:)),2);
        AdPolar_1{k}(:,i)=sum(AdPolar(:,Polar1{k}(i,:)),2);
    end
    [~,r_p1(k)]=size(Polar_1{k});
    for j=1:1:r_p1(k)
        maxbind_po{k}(j)=max(Polar_1{k}(:,j));
        minbind_po{k}(j)=min(Polar_1{k}(:,j));
        po{k}(j)=(maxbind_po{k}(j)-minbind_po{k}(j))+1;
    end        
end

for x=1:1:9
    q2=1;
    opt_0=[];
    for i=1:1:r_p(x)
        for a=minbind_po{x}(i):1:maxbind_po{x}(i)
             flags=Polar_1{x}(:,i)<a;
            for f=1:50
                Fl=flags*1;
                Fl(clim(:,f),:)=[];
                Fl_1(:,f)=Fl;
            end
            averghc=sum(flags(1:(hc),:),1)/(hc)*100;
            averglc=(sum(flags((hc+1):end,:),1)/(lc)*100);
            stdhc=std(sum(Fl_1(1:(hc-nhc)),1)/(hc-nhc)*100);
            stdlc=std(sum(Fl_1((hc-nhc+1):end),1)/(lc-nlc)*100);
            if averghc-stdhc<high && averglc+stdlc>low
               mAb=sum(flags)/(hc+lc)*100;
               highspec=averghc;
               lowspec=averglc;
               flag=a;
               Flags=flags';
               table(1,1)=sum(flags(1:hc)<1);
               table(2,1)=hc-table(1,1);
               table(1,2)=sum(flags(hc+1:end)<1);
               table(2,2)=lc-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=hc+lc;
               pvalue=TwoTailed(table);
               
               Adflags=AdPolar_1{x}(:,i)<a;
               Pflags=PPolar_1{x}(:,i)<a;
               Aflags=APolar_1{x}(:,i)<a;
               Cflags=CPolar_1{x}(:,i)<a;
               Eflags=EPolar_1{x}(:,i)<a;
               Bflags=BPolar_1{x}(:,i)<a;
               
               table(1,1)=sum(Adflags(1:ha)<1);
               table(2,1)=ha-table(1,1);
               table(1,2)=sum(Adflags(ha+1:end)<1);
               table(2,2)=la-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=ha+la;
               pvalue1=TwoTailed(table);
               
               for f=1:50
                   Ad=Adflags;
                   Pflags1=Pflags;
                   Aflags1=Aflags;
                   Cflags1=Cflags;
                   Eflags1=Eflags;
                   Bflags1=Bflags;
                   
                   Ad(adi(:,f),:)=[];
                   Pflags1(clis_p(:,f),:)=[];
                   Aflags1(clis_a(:,f),:)=[];
                   Cflags1(clis_c(:,f),:)=[];
                   Eflags1(clis_e(:,f),:)=[];
                   Bflags1(clis_b(:,f),:)=[];
                   
                   Ad_1(:,f)=Ad;
                   Pflags1_1(:,f)=Pflags1;
                   Aflags1_1(:,f)=Aflags1;
                   Cflags1_1(:,f)=Cflags1;
                   Eflags1_1(:,f)=Eflags1;
                   Bflags1_1(:,f)=Bflags1;
               end
               
               AdflagH=sum(Ad_1(1:(ha-nha),:),1)/(ha-nha)*100;
               AdflagL=sum(Ad_1((ha-nha+1):end,:),1)/(la-nla)*100;
               
               Pflag1H=sum(Pflags1_1(1:(gp-ngp),:),1)/(gp-ngp)*100;
               Aflag1H=sum(Aflags1_1(1:(ga-nga),:),1)/(ga-nga)*100;
               Cflag1H=sum(Cflags1_1(1:(gc-ngc),:),1)/(gc-ngc)*100;
               Eflag1H=sum(Eflags1_1(1:(ge-nge),:),1)/(ge-nge)*100;
               Bflag1H=sum(Bflags1_1(1:(gb-ngb),:),1)/(gb-ngb)*100;
               
               Pflag1L=sum(Pflags1_1((gp-ngp+1):end,:),1)/(bp-nbp)*100;
               Aflag1L=sum(Aflags1_1((ga-nga+1):end,:),1)/(ba-nba)*100;
               Cflag1L=sum(Cflags1_1((gc-ngc+1):end,:),1)/(bc-nbc)*100;
               Eflag1L=sum(Eflags1_1((ge-nge+1):end,:),1)/(be-nbe)*100;
               Bflag1L=sum(Bflags1_1((gb-ngb+1):end,:),1)/(bb-nbb)*100;
               
               avergha=sum(Adflags(1:ha,:),1)/ha*100;
               avergla=sum(Adflags(ha+1:end,:),1)/la*100;
               averghcp=sum(Pflags(1:gp,:),1)/(gp)*100;
               averghca=sum(Aflags(1:ga,:),1)/(ga)*100;
               averghcc=sum(Cflags(1:gc,:),1)/(gc)*100;
               averghce=sum(Eflags(1:ge,:),1)/(ge)*100;
               averghcb=sum(Bflags(1:gb,:),1)/(gb)*100;
               averglcp=sum(Pflags(gp+1:end,:),1)/(bp)*100;
               averglca=sum(Aflags(ga+1:end,:),1)/(ba)*100;
               averglcc=sum(Cflags(gc+1:end,:),1)/(bc)*100;
               averglce=sum(Eflags(ge+1:end,:),1)/(be)*100;
               averglcb=sum(Bflags(gb+1:end,:),1)/(bb)*100;
              
               stdha=std(AdflagH);
               stdla=std(AdflagL);
               stdhcp=std(Pflag1H);
               stdlcp=std(Pflag1L);
               stdhca=std(Aflag1H);
               stdlca=std(Aflag1L);
               stdhcc=std(Cflag1H);
               stdlcc=std(Cflag1L);
               stdhce=std(Eflag1H);
               stdlce=std(Eflag1L);
               stdhcb=std(Bflag1H);
               stdlcb=std(Bflag1L);  
               if pvalue<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10 
                  opt{x+13}(q2,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
                  opt_0=zeros(1,13-x);
                  opt_1{x+13}(q2,:)=[100,opt_0,Polar_num(Polar1{x}(i,:))];
                  opt_2{x+13}(q2,:)=Flags;
                  q2=q2+1;
               elseif pvalue1<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10  
                  opt{x+13}(q2,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
                  opt_0=zeros(1,13-x);
                  opt_1{x+13}(q2,:)=[100,opt_0,Polar_num(Polar1{x}(i,:))];
                  opt_2{x+13}(q2,:)=Flags;
                  q2=q2+1;
               end
            end
        end
    end
end

%% single amino acid
Single=netcharge;
PSingle=Pnetcharge;
ASingle=Anetcharge;
CSingle=Cnetcharge;
ESingle=Enetcharge;
BSingle=Bnetcharge;
AdSingle=ABnetcharge;
y=1;
for a=0:1:max_nc
    flags=Single>a;
    for f=1:50
        Fl=flags*1;
        Fl(clim(:,f),:)=[];
        Fl_1(:,f)=Fl;
    end
    averghc=(sum(flags(1:(hc),:),1)/(hc)*100);
            averglc=(sum(flags((hc+1):end,:),1)/(lc)*100);
            stdhc=std(sum(Fl_1(1:(hc-nhc)),1)/(hc-nhc)*100);
            stdlc=std(sum(Fl_1((hc-nhc+1):end),1)/(lc-nlc)*100);
            if averghc-stdhc<high && averglc+stdlc>low
               mAb=sum(flags)/(hc+lc)*100;
               highspec=averghc;
               lowspec=averglc;
               flag=a;
               Flags=flags';
               table(1,1)=sum(flags(1:hc)<1);
               table(2,1)=hc-table(1,1);
               table(1,2)=sum(flags(hc+1:end)<1);
               table(2,2)=lc-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=hc+lc;
               pvalue=TwoTailed(table);
       
       Adflags=AdSingle>a;
       Pflags=PSingle>a;
       Aflags=ASingle>a;
       Cflags=CSingle>a;
       Eflags=ESingle>a;
       Bflags=BSingle>a;
       
       table(1,1)=sum(Adflags(1:ha)<1);
               table(2,1)=ha-table(1,1);
               table(1,2)=sum(Adflags(ha+1:end)<1);
               table(2,2)=la-table(1,2);
               table(3,1)=table(1,1)+table(2,1);
               table(3,2)=table(1,2)+table(2,2);
               table(1,3)=table(1,1)+table(1,2);
               table(2,3)=table(2,1)+table(2,2);
               table(3,3)=ha+la;
               pvalue1=TwoTailed(table);
       for f=1:50
                   Ad=Adflags;
                   Pflags1=Pflags;
                   Aflags1=Aflags;
                   Cflags1=Cflags;
                   Eflags1=Eflags;
                   Bflags1=Bflags;
                   
                   Ad(adi(:,f),:)=[];
                   Pflags1(clis_p(:,f),:)=[];
                   Aflags1(clis_a(:,f),:)=[];
                   Cflags1(clis_c(:,f),:)=[];
                   Eflags1(clis_e(:,f),:)=[];
                   Bflags1(clis_b(:,f),:)=[];
                   
                   Ad_1(:,f)=Ad;
                   Pflags1_1(:,f)=Pflags1;
                   Aflags1_1(:,f)=Aflags1;
                   Cflags1_1(:,f)=Cflags1;
                   Eflags1_1(:,f)=Eflags1;
                   Bflags1_1(:,f)=Bflags1;
               end
               
               AdflagH=sum(Ad_1(1:(ha-nha),:),1)/(ha-nha)*100;
               AdflagL=sum(Ad_1((ha-nha+1):end,:),1)/(la-nla)*100;
               
               Pflag1H=sum(Pflags1_1(1:(gp-ngp),:),1)/(gp-ngp)*100;
               Aflag1H=sum(Aflags1_1(1:(ga-nga),:),1)/(ga-nga)*100;
               Cflag1H=sum(Cflags1_1(1:(gc-ngc),:),1)/(gc-ngc)*100;
               Eflag1H=sum(Eflags1_1(1:(ge-nge),:),1)/(ge-nge)*100;
               Bflag1H=sum(Bflags1_1(1:(gb-ngb),:),1)/(gb-ngb)*100;
               
               Pflag1L=sum(Pflags1_1((gp-ngp+1):end,:),1)/(bp-nbp)*100;
               Aflag1L=sum(Aflags1_1((ga-nga+1):end,:),1)/(ba-nba)*100;
               Cflag1L=sum(Cflags1_1((gc-ngc+1):end,:),1)/(bc-nbc)*100;
               Eflag1L=sum(Eflags1_1((ge-nge+1):end,:),1)/(be-nbe)*100;
               Bflag1L=sum(Bflags1_1((gb-ngb+1):end,:),1)/(bb-nbb)*100;
               
               avergha=sum(Adflags(1:ha,:),1)/ha*100;
               avergla=sum(Adflags(ha+1:end,:),1)/la*100;
               averghcp=sum(Pflags(1:gp,:),1)/(gp)*100;
               averghca=sum(Aflags(1:ga,:),1)/(ga)*100;
               averghcc=sum(Cflags(1:gc,:),1)/(gc)*100;
               averghce=sum(Eflags(1:ge,:),1)/(ge)*100;
               averghcb=sum(Bflags(1:gb,:),1)/(gb)*100;
               averglcp=sum(Pflags(gp+1:end,:),1)/(bp)*100;
               averglca=sum(Aflags(ga+1:end,:),1)/(ba)*100;
               averglcc=sum(Cflags(gc+1:end,:),1)/(bc)*100;
               averglce=sum(Eflags(ge+1:end,:),1)/(be)*100;
               averglcb=sum(Bflags(gb+1:end,:),1)/(bb)*100;
              
               stdha=std(AdflagH);
               stdla=std(AdflagL);
               stdhcp=std(Pflag1H);
               stdlcp=std(Pflag1L);
               stdhca=std(Aflag1H);
               stdlca=std(Aflag1L);
               stdhcc=std(Cflag1H);
               stdlcc=std(Cflag1L);
               stdhce=std(Eflag1H);
               stdlce=std(Eflag1L);
               stdhcb=std(Bflag1H);
               stdlcb=std(Bflag1L);
       if pvalue<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10  
          opt{45}(y,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
          opt_1{45}(y,:)=[710,0,0,0,0,0,0,0,0,0,0,0,0,1];
          opt_2{45}(y,:)=[Flags];
          y=y+1;
       elseif pvalue1<0.05 && averghcp+stdhcp>=avergha-stdha && averglcp-stdlcp<=avergla+stdla &&  avergha<avergla-10 && averghcp<averglcp-10 && averghca<averglca-10 && averghcc<averglcc-10 && averghce<averglce-10 && averghcb<averglcb-10  
          opt{45}(y,:)=[flag,pvalue,mAb,100-averghc,100-averglc,100-averghcp,100-averglcp,100-averghca,100-averglca,100-averghcc,100-averglcc,100-averghce,100-averglce,100-averghcb,100-averglcb,100-avergha,100-avergla,stdhc,stdlc,stdhcp,stdlcp,stdhca,stdlca,stdhcc,stdlcc,stdhce,stdlce,stdhcb,stdlcb,stdha,stdla,sum(flags),sum(Pflags),sum(Aflags),sum(Cflags),sum(Eflags),sum(Bflags),pvalue1];
          opt_1{45}(y,:)=[710,0,0,0,0,0,0,0,0,0,0,0,0,1];
          opt_2{45}(y,:)=[Flags];
          y=y+1;
       end
    end
end
               
exist opt               
if ans~=0
   [~,r]=size(opt);
   Res=[];
   for i=1:r
       Res=[Res;opt_1{i} opt{i} opt_2{i}];
   end
   Res(isnan(Res(:,1)),:) = [] ;

   [e,~]=size(Res);
   Fx=Res(:,1:14);
   Px=Res(:,15:19);
   q=0;
   for i=1:e-1
       for j=i+1:e
           p=ismember(Fx(i,:),Fx(j,:));
           p1=sum(p);
           if p1==14 && Px(i,1)==Px(j,1)&& Px(i,2)==Px(j,2) && Px(i,3)==Px(j,3)
              q=q+1;
              Row(q)=j;
           end
       end
   end

   if q~=0
      Row=unique(Row);
      Res(Row,:) = []; 
      Res(:,18)=-1*Res(:,18);
      Res(:,20)=-1*Res(:,20);
      Res(:,22)=-1*Res(:,22);
      Res(:,24)=-1*Res(:,24);
      Res(:,26)=-1*Res(:,26);
      Res(:,28)=-1*Res(:,28);
      Res(:,30)=-1*Res(:,30);
      Res=sortrows(Res,[16,18,30,20,22,24,26,28]);
      Res(:,18)=-1*Res(:,18);
      Res(:,20)=-1*Res(:,20);
      Res(:,22)=-1*Res(:,22);
      Res(:,24)=-1*Res(:,24);
      Res(:,26)=-1*Res(:,26);
      Res(:,28)=-1*Res(:,28);
      Res(:,30)=-1*Res(:,30);
      [res,~]=size(Res);
      Sig=zeros(res,1);
      for i=1:res
          pc=i/res*0.05;
          if Res(i,16)<pc
             Sig(i,:)=1;
          end
      end
      s=find(Sig==1);
      [rs,~]=size(s);
      if rs==0
         disp('No results')
      else
      max1=max(s);
      Res_1=Res(1:max1,:); % results from FDR
      end
   else
      Res(:,18)=-1*Res(:,18);
      Res(:,20)=-1*Res(:,20);
      Res(:,22)=-1*Res(:,22);
      Res(:,24)=-1*Res(:,24);
      Res(:,26)=-1*Res(:,26);
      Res(:,28)=-1*Res(:,28);
      Res(:,30)=-1*Res(:,30);
      Res=sortrows(Res,[16,18,30,20,22,24,26,28]);
      Res(:,18)=-1*Res(:,18);
      Res(:,20)=-1*Res(:,20);
      Res(:,22)=-1*Res(:,22);
      Res(:,24)=-1*Res(:,24);
      Res(:,26)=-1*Res(:,26);
      Res(:,28)=-1*Res(:,28);
      Res(:,30)=-1*Res(:,30);
      [res,~]=size(Res);
      Sig=zeros(res,1);
      for i=1:res
          pc=i/res*0.05;
          if Res(i,16)<pc
          Sig(i,:)=1;
          end
      end
      s=find(Sig==1);
      [rs,~]=size(s);
      if rs==0
         disp('No results')
      else
         max1=max(s);
         Res_1=Res(1:max1,:); % results from FDR
      end
   end
   exist Res_1
   if ans ~=0 
      [res1,~]=size(Res_1);
      if res1~=0
         Sign_v=Res_1(:,1);
         Flag_AA=Res_1(:,2:14);
         other_info=Res_1(:,15:end);
         for i=1:res1
             F=nonzeros(Flag_AA(i,:))';
             if Sign_v(i)==510
                Hydro_AA='GAVLIMFWPYRKH';
                Flag_AA1{i}=Hydro_AA(F);
                Sign1{i}='>';
             elseif Sign_v(i)==100
                Pola_AA='GASTQNYDE';
                Flag_AA1{i}=Pola_AA(F);
                Sign1{i}='<';
             elseif Sign_v(i)==710
                Flag_AA1{i}='charge';
                Sign1{i}='>';
             end
         end
      end
      Flag=cell2table(Flag_AA1');
      Sign=cell2table(Sign1');
      Flag.Properties.VariableNames{'Var1'}='Flag';
      Sign.Properties.VariableNames{'Var1'}='Sign';
      info=array2table(other_info);
      R=cell(res1,1);
      R(:)={Region(re2)};
      Singleflag=[R,Flag,Sign,info];
   else
      Singleflag=[];
   end
   
  else
   Singleflag=[];
end
re2
  Finalflag=[Finalflag;Singleflag];
end


