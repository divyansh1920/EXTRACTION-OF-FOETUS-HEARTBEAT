% clc
% clear all;
% close all;
load foetal_ecg.dat %%%loading the input signal%%%
S= foetal_ecg; % source signal
Fs=500; % sampling Frequency
t=S(:,1); % Time samples
%%%PLOTTING INPUT SIGNALS%%%%%%
figure
d1=S(:,2); %%%Abdominal signal 1
subplot(3,1,1);
plot(t,d1,'r');
xlabel('time period(sec)');
ylabel('amplitude(mV)');
title('Abdominal signal 1');
d2=S(:,3); %%%Abdominal signal 2
subplot(3,1,2);
plot(t,d2,'r');
xlabel('time period(sec)');
ylabel('amplitude(mV)');
title('Abdominal signal 2');
d3=S(:,4); %%%Abdominal signal 3
subplot(3,1,3);
plot(t,d3,'r');
xlabel('time period(sec)');
ylabel('amplitude(mV)');
title('Abdominal signal 3');
d4=S(:,5); %%%Abdominal signal 4
figure
subplot(2,1,1);
plot(t,d4,'r');
xlabel('time period(sec)');
ylabel('amplitude(mV)');
title('Abdominal signal 4');
d5=S(:,6); %%%Abdominal signal 5
subplot(2,1,2);
plot(t,d5,'r');
xlabel('time period(sec)');
ylabel('amplitude(mV)');
title('Abdominal signal 5');
figure
x1=S(:,7); %Thoracic signal 1
subplot(3,1,1);
plot(t,x1,'b');
xlabel('time period(sec)');
ylabel('amplitude(mV)');
title('Thoracic signal 1');
x2=S(:,8); %Thoracic signal 2
subplot(3,1,2);
plot(t,x2,'b');
xlabel('time period(sec)');
ylabel('amplitude(mV)');
title('Thoracic signal 2');
x3=S(:,9); %Thoracic signal 3
subplot(3,1,3);
plot(t,x3,'b');
xlabel('time period(sec)');
ylabel('amplitude(mV)');
title('Thoracic signal 3');
d=(d1+d2+d3+d4+d5)/5; %%% AVERAGE OF ABDOMINAL SIGNALS
x=(x1+x2+x3)/3; %%% AVERAGE OF THORACIC SIGNALS
figure
subplot(2,1,1);
plot(t,d,'r');
xlabel('time period(sec)');
ylabel('amplitude(mV)');
title('Average of abdominal signals[Mother+Fetus]');
subplot(2,1,2);
plot(t,x,'b');
xlabel('time period(sec)');
ylabel('amplitude(mV)');
title('Average of Thoracic signals(Mother)');

%SISO Implementation
%%% Generating ANC using LMS Algorithm
p=12;%order of filter
mu=0.0000005; % Step size
[A1,L,yl]=lms(x,d,mu,p);%calling LMS function


%%% Generating ANC using NLMS Algorithm
beta=0.001;%normalized step size
p=12;%order of filter
[A,LN,yn]=nlms(x,d,beta,p);%calling NLMS function


% %%% Generating ANC using LLMS Algorithm
gama= 0.001;%leakage coefficient
p=12;%order of filter
mu=0.0000009; % Step size
[AL,LL,yll]= llms(x,d,mu,gama,p);%calling LLMS function

%%%%plotting of lms signal%%%%
figure
subplot(2,1,1)
plot(t,L,'r');
legend('SISO-LMS');
title('Plot of the error signal-SISO');
ylabel('amplitude(mV)');
xlabel('time(SEC)-->');
axis([0 5 -30 30]);
%%%%%plotting of filtered signal of lms in siso%%%
subplot(2,1,2)
plot(t,yl,'r');
legend('SISO-LMS');
title('Plot of the filtered output-SISO');
ylabel('amplitucde(mV)');
xlabel('time(SEC)-->');
axis([0 2 -10 10]);

%plotting nlms signals
figure
subplot(2,1,1)
plot(t,LN,'r');
legend('SISO-NLMS');
title('Plot of the error signal-SISO');
ylabel('amplitude(mV)');
xlabel('time(SEC)-->');
axis([0 5 -30 30]);
%%%%%plotting of filtered signal of nlms in siso%%%
subplot(2,1,2)
plot(t,yn,'r');
legend('SISO-NLMS');
title('Plot of the filter output-SISO');
ylabel('amplitucde(mV)');
xlabel('time(SEC)-->');
axis([0 5 -30 30]);

%%%%plotting of llms signal%%%%
figure
subplot(2,1,1)
plot(t,LL,'r');
legend('SISO-LLMS');
title('Plot of the error signal-SISO');
ylabel('amplitude(mV)');
xlabel('time(SEC)-->');
axis([0 5 -30 30]);
%%%%%plotting of filtered signal of lms in siso%%%
subplot(2,1,2)
plot(t,yll,'r');
legend('SISO-LLMS');
title('Plot of the filter output-SISO');
ylabel('amplitucde(mV)');
xlabel('time(SEC)-->');
axis([0 2 -10 10]);

%COMPARISION
%for error
figure
plot(t,L,'r-',t,LN,'b-',t,LL,'g--');
legend('LMS-error','NLMS-error','LLMS-error');
title('Plot of the LMS,NLMS,LLMS output[Foetus signal]-error-SISO ');
ylabel('amplitude(mV)');
xlabel('time(SEC)-->');
axis([0 5 -30 30]);

%for output
figure
plot(t,yl,'r-',t,yn,'b-',t,yll,'g--');
legend('LMS-output','NLMS-output','LLMS-output');
title('Plot of the LMS,NLMS,LLMS output -SISO');
ylabel('amplitude(mV)');
xlabel('time(SEC)-->');
axis([0 5 -30 30]);


%% ERROR BETWEEN SISO and MISO%%
figure
% % LMS ERROR:
yl=reshape(yl,1,2500);
lms_error = abs(ym1-yl);
subplot(311)
plot(lms_error);
% NLMS ERROR:
nlms_error = abs(ym2-yn);
subplot(312)
plot(nlms_error);
% % LLMS ERROR:
yll=reshape(yll,1,2500);
llms_error = abs(ym3-yll);
subplot(313)
plot(llms_error);

%% SNR MISO%%
% % LMS ERROR:
lms_snr_miso = 20 * log10(rms(ym1)/rms(d) );

% NLMS ERROR:
nlms_snr_miso = 20 * log10(rms(ym2)/rms(d) );

% % LLMS ERROR:
llms_snr_miso = 20 * log10(rms(ym3)/rms(d) );

%% SNR SISO%%
% % LMS ERROR:
lms_snr_siso = 20 * log10(rms(yl)/rms(d) );

% NLMS ERROR:
nlms_snr_siso = 20 * log10(rms(yn)/rms(d) );

% % LLMS ERROR:
llms_snr_siso = 20 * log10(rms(yll)/rms(d) );
%%
%LMS Function
% function [w,y,e,W] = lms(x,d,mu_step,M)
% N = length(x); % number of data samples
% y = zeros(N,1); % initialize filter output vector
% w = zeros(M,1); % initialize filter coefficient vector
% e = zeros(N,1); % initialize error vector
% W = zeros(M,N); % filter coefficient matrix for coeff. history
% for n = 1:N
%   if n <= M % assume zero-samples for delayed data that isn't available
%       k = n:-1:1;
%       x1 = [x(k); zeros(M-numel(k),1)];
%   else
%       x1 = x(n:-1:n-M+1); % M samples of x in reverse order
%   end
%   y(n) = w'*x1; % filter output
%   e(n) = d(n) - y(n); % error
%   w = w + mu_step*e(n)'*x1; % update filter coefficients
%   W(:,n) = w; % store current filter coefficients in matrix
% end
% end
function [A,E,Y] = lms(x,d,mu,nord)
X=convm(x,nord);
[M,N]=size(X);
if nargin < 5, a0 = zeros(1,N); end
a0=a0(:).';
Y(1)=a0*X(1,:).';
E(1)=d(1) - Y(1);
A(1,:) = a0 + mu*E(1)*conj(X(1,:));
if M>1
for k=2:M-nord+1;
Y(k,:)=A(k-1,:)*X(k,:).';%ouput equation
E(k,:) = d(k) - Y(k,:);%error signal
A(k,:)=A(k-1,:)+mu*E(k)*conj(X(k,:));%update equation
end
end
end

% NLMS Function:
%%NLMS CALGORITHM FOR THE SOURSE CODE%%
function [A,E,Y] = nlms(x,d,beta,nord)
X=convm(x,nord);
[M,N]=size(X);
if nargin < 5, a0 = zeros(1,N); end%initialization
a0=a0(:).';
Y(1)=a0*X(1,:).';
E(1)=d(1) - a0*X(1,:).';
DEN=X(1,:)*X(1,:)'+0.0001;
A(1,:) = a0 + beta/DEN*E(1)*conj(X(1,:));
    if M>1
        for k=2:M-nord+1;
            Y(k)=A(k-1,:)*X(k,:).';%output equation
            E(k) = d(k) - A(k-1,:)*X(k,:).';%error signal
            DEN=X(k,:)*X(k,:)'+0.0001;%normalizing the input signal
            A(k,:)=A(k-1,:)+ beta/DEN*E(k)*conj(X(k,:));%update equation
        end
    end
end

%LLMS Function:
%%LLMS FUNCTION OF THE SOURCE CODE%%
function [A,E,Y]= llms(x,d,mu,gama,nord,a0)
X=convm(x,nord);

[M,N]=size(X);
if nargin < 6, a0 = zeros(1,N); end
a0=a0(:).';
Y(1)=a0*X(1,:).';
E(1)=d(1) - Y(1);
A(1,:)=(1-mu*gama)*a0+mu*E(1)*conj(X(1,:));
if M>1
for k=2:M-nord+1
Y(k,:)=A(k-1,:)*X(k,:).';%output signal
E(k,:) = d(k) - Y(k,:);%error signal
A(k,:)=(1-mu*gama)*A(k-1,:)+mu*E(k)*conj(X(k,:));%update eqtn
end
end
end



