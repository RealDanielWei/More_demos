function y=min_tstep(z)
X=sqrt(real(z).*real(z)-3*imag(z).*imag(z));
Z=abs(z);
y=sqrt((2*real(z)-2*X)./(3*Z.*Z));
end