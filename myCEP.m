function [cep_coef,bank_coef,freq_coef,freq_energy] = myCEP(data,fs,...
    window,frame,nfft,banks,conds,limL,limH)

duration = length(data);
[x,mc,na,nb] = melbankm(banks,nfft,fs,limL,limH,conds);
elements = floor(duration/frame)-1;
% freq coefficients
freq_coef = zeros(nfft/2+1,elements);
% bank coefficients
bank_coef = zeros(banks,elements);
% cepstrum coefficents
cep_coef = zeros(banks,elements);
% E_f!
freq_energy = zeros(1,elements);
index = 1:frame:duration-window;
for i=1:numel(index)
    f=rfft(data(index(i):index(i)+window-1),nfft);
    b=log(x*(f(na:nb).*conj(f(na:nb)))');
    c = dct(b);
    % store values if required
    freq_coef(:,i) = f;
    bank_coef(:,i) = b;
    cep_coef(:,i) = c;
    freq_energy(i) = log(sum((abs(b(2:end)).^2)));
end

end