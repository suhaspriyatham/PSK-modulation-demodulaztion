clear all;
clc;

%PSK with i/p arguments as l:length of sequence, 8 symbols, Es:Energy of
%symbol, SNR from 0dB to 20dB 
% AWGN channel

l=3000000;
k=randi([1 8],1,l);
figure;
histogram(k);
x=zeros(2,l);
Es=1;
SNR=[0:20];

%8 symbols
for i=1:8
    symbols(:,i)=[Es*cos(2*pi*i/8); Es*sin(2*pi*i/8)];
end

%sequence of symbols
for i=1:l
    Sm(:,i)=[Es*cos(2*pi*k(i)/8); Es*sin(2*pi*k(i)/8)];
end


%setting the SNR and variance(sigma^2)
for snrloop=1:length(SNR)
    
    sigma=sqrt((Es/6)*10^(-SNR(snrloop)/10));
    No(snrloop)=2*sigma^2;
    Eb=Es/3;
    snr(snrloop)=Eb/No(snrloop);
    %generating white gaussian noise
    mean=0;
    n=normrnd(mean,sigma,[2,l]);
%     figure;
%     histogram(n(1,:));
%     figure;
%     histogram(n(2,:));


    % received vector through AWGN channel
    R=Sm+n;
    figure;
    plot(R(1,:),R(2,:),'.');
    hold on;
    plot(symbols(1,:),symbols(2,:),'r*');
    xlabel('In phase Component');
    ylabel('Quadrature Component');
    title(['PSK constellation for SNR ',num2str(SNR(snrloop)),'dB']);
    %axis([-4 4 -4 4]);
    grid on;
    hold off;

    
    %Maximum likelihood
    for i=1:l
        for j=1:8
            d(j)=(symbols(1,j)-R(1,i))^2 + (symbols(2,j)- R(2,i))^2;
        end
        [dmin,S]=min(d);
        M(i)=S;
    end

    %comparing sent sequence k and received sequence M
    if M==k
        X=['SNR ',num2str(SNR(snrloop)),'dB No error'];
        disp(X);
    else
        X=['SNR ',num2str(SNR(snrloop)),'dB Received error sequence due to noise'];
        disp(X);
    end

    cr=0;
    w=0;
    for i=1:l
        if M(i)==k(i)
            cr=cr+1;
        else
            w=w+1;
        end
    end 
    SER(snrloop)=w/l;

end

figure;
semilogy(SNR,SER);
xlabel('SNR in dB');
ylabel('Symbol error rate');
title('SNR vs SER');
hold on;

axis([0 25 0 1]);
% figure;
% plot(SNR,SER);
% xlabel('SNR in dB');
% ylabel('Symbol error rate');
% title('SNR vs SER');
% figure;
% loglog(snr,SER);
% axis([0 25 10^-7 1]);

%theoritical values solving
Dmin=2*sin((360/16)*pi/180)
Nmin=2;
m=8;    %number of symbols
Pemax=7*.5*erfc(sqrt(3*10.^(SNR/10)*(1-cos(2*pi/m))/2));%(m-1)*qfunc(sqrt((Dmin^2)/2*No));
Pemin=.5*erfc(sqrt(3*10.^(SNR/10)*(1-cos(2*pi/m))/2));%(Nmin/m)*qfunc(sqrt((Dmin^2)/2*No));

semilogy(SNR,Pemin,'r'); 
semilogy(SNR,Pemax,'g');
hold off;