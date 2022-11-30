%Unfortunately, saving an octave file with Greek characters in the comments results in the characters
%turning into gibberish. As a result, the original comments have been lost.

%I am not sure what this script was for, as the original report was lost. Judging from the context, and
%the few things I remember, it was probably a script used in Introduction to Chemical Engineering (EXM)
%where we had to calculate the moisture in a product.
%The price of the product was different depending on the moisture content
%I have attempted to recover as much as possible.

clear;clc;
%Exercise 3
A=[0.1132 0.0453 0.0404 0.0379 0.0451 0.0459 0.0405 0.0372 0.0387 0.0402];
%B was the moisture content at the samples, determined randomly
B=rand(1,100);
B=B.*0.1; 
MoistureMean=mean(A); %The array A is used to determine the mean moisture in the samples.
disp('MoistureMean is: '), disp(MoistureMean);
C=MoistureMean*1.1<B; %The goal is to determine the days where the moisture content exceeds the average by
% 10 percent. C is a vector which details the days where the moisture content is greater than 10% the mean
d=sum(C); %The result of the above vector is 1 and zero, representing a "True/False" state.
%The variable d represents the days where the moisture is greater than the average plus 10%
disp('The amount of days where the moisture is larger than the mean by 10% is: ')
disp(d)
%Exercise 4
%Calculation of the gain from selling the products
moistprof=d*(220*10^3)*(0.11); %Gain from the off-spec product
driedprof=(100-d)*(200*10^3)*(0.18); %Gain from the product
%Displying to the user
disp('The profit from the mixture with moisture above the limit is (in euros):');
disp(moistprof);
disp('The profit from the mixture with moisture below the limit is (in euros):');
disp(driedprof);


