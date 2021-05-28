clear all; close all; clc;

%test 1
%Correct values
%From the values when data is defined.
cw = [ 1     0     0     0;
     0     1     0     0;
     0     0     1     0];
cwH = MakeHomoPts(cw);
tFromCtoW = [
    0.865815567821920,0.459788742279509,0.197377088311267,0.593610836073900;
    -0.498738747594092,0.761262407411635,0.414414295977007,0.281518366646709;
    0.0402872705141819,-0.457245950812984,0.888427305016931,6.12624262406678;
    0,0,0,1];
ccH = (tFromCtoW * MakeHomoPts(cw));
cc = ccH(1:3,:);
TheWH = [   
    -1.6830    0.0254   -1.6426    1.2273   -0.1736   -0.2507   -0.2042    0.8618    0.8189    1.0206;
   -0.7113    1.2391   -1.1975   -0.4020    0.6236    0.5754    1.3177   -1.4390    0.8344   -0.8406;
   -0.2195   -0.2226    2.4573    0.1638    0.0193   -2.0083    0.3477   -1.0972    0.2345    0.3251;
    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000];
TheCH = [
       -1.2339    1.1414   -0.8941    1.5037    0.7339    0.2448    1.0913    0.4616    1.7326    1.1549;
    0.4884    1.1199    1.2075   -0.5687    0.8508    0.0123    1.5306   -1.6985    0.6055   -0.7327;
    6.1886    5.3630    8.7907    6.5050    5.8513    4.0688    5.8244    5.8441    5.9860    6.8405;
    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000];
TheW = TheWH(1:3,:);
TheC = TheCH(1:3,:);
%This is from TheC,TheW
[TcToW,aIsBiggerBy] = GetScaleTransFromA2B(TheC,TheW);
errorInTcToW = ErrorIn(tFromCtoW,TcToW)    

%Solving using our method
eigVector = [-0.114375735144126;0.0170236338693208;-0.483273009244033;-0.0825552770958948;-0.0817230778535370;-0.444281164368532;-0.0619899880529638;-0.0545404754127483;-0.549742027272243;-0.0465214796609366;-0.0220626885027153;-0.480115682386221];
beta = 12.759930260180505;
ccour = -reshapeVector(beta*eigVector)';
unscaledCC = reshapeVector(eigVector)';
errorInCCourCC = ErrorIn(ccour,cc)

xcour = [
    -1.2339    1.1414   -0.8941    1.5037    0.7339    0.2448    1.0913    0.4616    1.7326    1.1549;
    0.4884    1.1199    1.2075   -0.5687    0.8508    0.0123    1.5306   -1.6985    0.6055   -0.7327;
    6.1886    5.3630    8.7907    6.5050    5.8513    4.0688    5.8244    5.8441    5.9860    6.8405
];

%Our solution using cc and cw
[tctow,aIsBiggerBy] = GetScaleTransFromA2B(ccour,cw)
errorInTcToW = ErrorIn(tctow,TcToW)
[tctowCC,aIsBiggerBy,ccScaled,cwScaled] = GetScaleTransFromA2B(unscaledCC,cw)
ErrorIn(tctowCC * cwH,MakeHomoPts(ccScaled))
ErrorIn(tctowCC * MakeHomoPts(cwScaled),MakeHomoPts(unscaledCC))
ErrorIn(tctowCC,TcToW)
ErrorIn(tctowCC * cwH,ccH)


% %Don't know what these are
% pointsA = randi([-10,10],4,50);
% pointsA(4,:)=1;
% % pointsA = [1 0 0 1; 0 1 0 1; 0 0 1 1; 0 0 0 1]'
% r = Rx(randi([-90,90])) * Ry(randi([-90,90])) * Rz(randi([-90,90]));
% % r = Rz(pi()/2);
% d = randi([-10,10],3,1);
% % d = [1 2 3 1]';
% T1 = MakeT(r,d,1);
% scale = randi([0,10]);
% % scale = 10;
% T = MakeT(r,d,scale);
% pointsBTRS = T1 * (scale * pointsA);
% pointsBTRS(4,:) = 1;
% pointsBSTR = scale*(T1*pointsA);
% pointsBSTR(4,:) = 1;
% pointsB = T * pointsA;
% 
% max(max(pointsBTRS-pointsBSTR))
% max(max(pointsBSTR-pointsB))
% 
% [Ti,scaleA2B] = invT(T)
% max(max(Ti * T - eye(4)))
% 
% [angB2A,axis,dist, scaleB2A] = Screws(T,1)
% [angA2B,axis,dist, scaleA2B] = Screws(Ti,1)
% 
% scaleB2A * scaleA2B - 1


function [error] = ErrorIn(a,b)
    error = max(max(abs(a-b)));
end