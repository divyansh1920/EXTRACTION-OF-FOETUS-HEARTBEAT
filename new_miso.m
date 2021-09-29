% MATLAB code for MISO system
clc;
clear all;
close all;
load('foetal_ecg.dat');
x=foetal_ecg;
% time signal;
timesig=x(:,1);
t=x(:,1);
% abdnomial signals
abdomin1=x(:,2);
abdomin2=x(:,3);
abdomin3=x(:,4);
abdomin4=x(:,5);
abdomin5=x(:,6);
%thoriad signals
thoirad1=x(:,7);
thoirad2=x(:,8);
thoirad3=x(:,9);
% figure
% subplot(5,1,1);
% plot(timesig,abdomin1);
% title('abdomin1');
% xlabel('time[s]');
% ylabel('amplitude mV');
% subplot(5,1,2);
% plot(timesig,abdomin2);
% title('abdomin2');
% ylabel('amplitude mV');
% xlabel('time');
% subplot(5,1,3);
% plot(timesig,abdomin3);
% title('abdomin3');
% xlabel('time');
% ylabel('amplitude mV');
% subplot(5,1,4);
% plot(timesig,abdomin4);
% title('abdomin4');
% xlabel('time');
% ylabel('amplitude mV');
% subplot(5,1,5);
% plot(timesig,abdomin5);
% title('abdomin5');
% xlabel('time');
% ylabel('amplitude mV');
% figure
% subplot(3,1,1);
% plot(timesig,thoirad1,'r');
% title('thoirad1');
% xlabel('time');
% ylabel('amplitude mV');
% subplot(3,1,2);
% plot(timesig,thoirad2,'r');
% title('thoirad2');
% xlabel('time');
% ylabel('amplitude mV');
% subplot(3,1,3);
% plot(timesig,thoirad3,'r');
% title('thoirad3');
% xlabel('time');
% ylabel('amplitude mV');
d=(abdomin1+abdomin2+abdomin3+abdomin4+abdomin5)/5;
a=thoirad1;
a1=thoirad2;
a2=thoirad3;
%% Applying for LMS Algorithm
mue= 0.00000002;
nord=12;
X=convm(a,nord);
X1=convm(a1,nord);
X2=convm(a2,nord);
%Applying LMS algorithm using lms basic function.
[A,A1,A2,E1,ym1] = lms1(X,X1,X2,d,mue,nord);
%% Applying for NLMS Algorithm
beta=0.005;
nord=12;
X=convm(a,nord);
X1=convm(a1,nord);
X2=convm(a2,nord);
%Applying NLMS algorithm using lms basic function.
[A,A1,A2,E2,ym2] = nlms1(X,X1,X2,d,beta,nord);
%% Applying for LLMS Algorithm
mu=0.0000002;
gammax=0.001;
nord=12;
X=convm(a,nord);
X1=convm(a1,nord);
X2=convm(a2,nord);
%Applying LMS algorithm using llms basic function.
[W,W1,W2,E3,ym3] = llms1(X,X1,X2,d,mu,gammax,nord);
%% Plotting signals
%%%%plotting of lms signal%%%%
figure
subplot(2,1,1)
plot(t,E1(1:2500),'r');
legend('MISO-LMS');
title('Plot of the error signal-MISO');
ylabel('amplitude(uVOLTS)');
xlabel('time(SEC)-->');
axis([0 5 -30 30]);
%%%%%plotting of filtered signal of lms in miso%%%
subplot(2,1,2)
plot(t,ym1(1:2500),'r');
legend('MISO-LMS');
title('Plot of the filtered output-MISO');
ylabel('amplitucde(uVOLTS)');
xlabel('time(SEC)-->');
axis([0 2 -10 10]);

%plotting nlms signals
figure
subplot(2,1,1)
plot(t,E2(1:2500),'r');
legend('MISO-NLMS');
title('Plot of the error signal-MISO');
ylabel('amplitude(uVOLTS)');
xlabel('time(SEC)-->');
axis([0 5 -30 30]);
%%%%%plotting of filtered signal of nlms in miso%%%
subplot(2,1,2)
plot(t,ym2(1:2500),'r');
legend('MISO-NLMS');
title('Plot of the filter output-MISO');
ylabel('amplitucde(uVOLTS)');
xlabel('time(SEC)-->');
axis([0 5 -30 30]);

%%%%plotting of llms signal%%%%
figure
subplot(2,1,1)
plot(t,E3(1:2500),'r');
legend('MISO-LLMS');
title('Plot of the error signal-MISO');
ylabel('amplitude(uVOLTS)');
xlabel('time(SEC)-->');
axis([0 5 -30 30]);
%%%%%plotting of filtered signal of lms in siso%%%
subplot(2,1,2)
plot(t,ym3(1:2500),'r');
legend('MISO-LLMS');
title('Plot of the filter output-MISO');
ylabel('amplitucde(uVOLTS)');
xlabel('time(SEC)-->');
axis([0 2 -10 10]);

%COMPARISION
%for error
figure
plot(t,E1(1:2500),'r-',t,E2(1:2500),'b-',t,E3(1:2500),'g--');
legend('LMS-error','NLMS-error','LLMS-error');
title('Plot of the MISO - LMS,NLMS,LLMS output[Foetus signal]-error ');
ylabel('amplitude(uVOLTS)');
xlabel('time(SEC)-->');
axis([0 5 -30 30]);

%for output
figure
plot(t,ym1(1:2500),'r-',t,ym2(1:2500),'b-',t,ym3(1:2500),'g--');
legend('LMS-output','NLMS-output','LLMS-output');
title('Plot of the LMS,NLMS,LLMS output -MISO');
ylabel('amplitude(uVOLTS)');
xlabel('time(SEC)-->');
axis([0 5 -30 30]);



% MISO Calling Functions:
% CONVM Function:
function[X] = convm(x,p)
N = length(x)+2*p-2; x = x(:);
xpad = [zeros(p-1,1);x;zeros(p-1,1)];
for i=1:p
X(:,i)=xpad(p-i+1 :N-i+1);
end;
end
% LMS Function:
function [A,A1,A2,E,y] = lms1(X,X1,X2,d,mu,nord,a0)
[M,N] = size(X);
[M1,N1] = size(X1);
[M2,N2] = size(X2);
if nargin < 7, a0 = zeros(1,N); end
a0 = a0(:).';
y1= zeros(1,M);
y2= zeros(1,M1);
y3= zeros(1,M2);
E=zeros(1,M);
E1=zeros(1,M1);
E2=zeros(1,M2);
A=zeros(size(X));
A1=zeros(size(X1));
A2=zeros(size(X2));
y1(1)= a0*X(1,:).';
y2(1)= a0*X1(1,:).';
y3(1)= a0*X2(1,:).';
E(1) = d(1) - a0*X(1,:).';
A(1,:) = a0 + mu*E(1)*conj(X(1,:));
A1(1,:) = a0 + mu*E(1)*conj(X1(1,:));
A2(1,:) = a0 + mu*E(1)*conj(X2(1,:));
if M>1
for k=2:M-nord+1;
y1(k) = A(k-1,:)*X(k,:).';
y2(k) = A1(k-1,:)*X1(k,:).';
y3(k) = A2(k-1,:)*X2(k,:).';
y(k)=(y1(k)+y2(k)+y3(k))/3;
E(k) = d(k) - y(k);
A(k,:) = A(k-1,:) + mu*E(k)*conj(X(k,:));
A1(k,:) = A1(k-1,:) + mu*E(k)*conj(X1(k,:));
A2(k,:) = A2(k-1,:) + mu*E(k)*conj(X2(k,:));
end;
end;
end
% NLMS Function:
function [A,A1,A2,E,y] = nlms1(X,X1,X2,d,beta,nord,a0)
[M,N] = size(X);
[M1,N1] = size(X1);
[M2,N2] = size(X2);
if nargin < 7, a0 = zeros(1,N); end
a0 = a0(:).';
y1= zeros(1,M);
y2= zeros(1,M1);
y3= zeros(1,M2);
E=zeros(1,M);
E1=zeros(1,M1);
E2=zeros(1,M2);
A=zeros(size(X));
A1=zeros(size(X1));
A2=zeros(size(X2));
y1(1)= a0*X(1,:).';
y2(1)= a0*X1(1,:).';
y3(1)= a0*X2(1,:).';
E(1) = d(1) - a0*X(1,:).';
DEN=X(1,:)*X(1,:)'+X1(1,:)*X1(1,:)' +X2(1,:)*X2(1,:)'+ 0.0001;
A(1,:) = a0 + beta/DEN*E(1)*conj(X(1,:));
A1(1,:) = a0 + beta/DEN*E(1)*conj(X1(1,:));
A2(1,:) = a0 + beta/DEN*E(1)*conj(X2(1,:));
if M>1
for k=2:M-nord+1;
y1(k) = A(k-1,:)*X(k,:).';
y2(k) = A1(k-1,:)*X1(k,:).';
y3(k) = A2(k-1,:)*X2(k,:).';
y(k)=(y1(k)+y2(k)+y3(k))/3;
E(k) = d(k) - y(k);
DEN=X(k,:)*X(k,:)'+X1(k,:)*X1(k,:)' +X2(k,:)*X2(k,:)'+ 0.0001;
A(k,:) = A(k-1,:) + beta/DEN*E(k)*conj(X(k,:));
A1(k,:) = A1(k-1,:) +beta/DEN*E(k)*conj(X1(k,:));
A2(k,:) = A2(k-1,:) + beta/DEN*E(k)*conj(X2(k,:));
end;
end;
end
% LLMS Function:
function [W,W1,W2,E,y] = llms1(X,X1,X2,d,mu,gammax,nord,a0)
[M,N] = size(X);
[M1,N1] = size(X1);
[M2,N2] = size(X2);
if nargin < 8, a0 = zeros(1,N);
end
a0 = a0(:).';
y1= zeros(1,M);
y2= zeros(1,M1);
y3= zeros(1,M2);
E=zeros(1,M);
E1=zeros(1,M1);
E2=zeros(1,M2);
W=zeros(size(X));
W1=zeros(size(X1));
W2=zeros(size(X2));
y1(1)= a0*X(1,:).';
y2(1)= a0*X1(1,:).';
y3(1)= a0*X2(1,:).';
E(1) = d(1) - a0*X(1,:).';
W(1,:) = (1 -mu*gammax)*a0 + mu*E(1)*conj(X(1,:));
W1(1,:) = (1 -mu*gammax)*a0 + mu*E(1)*conj(X1(1,:));
W2(1,:) = (1 -mu*gammax)*a0 + mu*E(1)*conj(X2(1,:));
if M>1
for k=2:M-nord+1
y1(k) = W(k-1,:)*X(k,:).';
y2(k) = W1(k-1,:)*X1(k,:).';
y3(k) = W2(k-1,:)*X2(k,:).';
y(k)=(y1(k)+y2(k)+y3(k))/3;
E(k) = d(k) - y(k);
W(k,:) = (1 -mu*gammax)*W(k-1,:) + mu*E(k)*conj(X(k,:));
W1(k,:) = (1 -mu*gammax)*W1(k-1,:) + mu*E(k)*conj(X1(k,:));
W2(k,:) = (1 -mu*gammax)*W2(k-1 ,:) + mu*E(k)*conj(X2(k,:));
end
end;
end