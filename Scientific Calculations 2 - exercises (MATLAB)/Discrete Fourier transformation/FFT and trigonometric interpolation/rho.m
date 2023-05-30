function o = rho(x,n)
  
  pom=dec2bin(x,n);
  obrnutistring=fliplr(pom);
  o=bin2dec(obrnutistring);
  
end