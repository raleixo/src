function string = tag2str(i)
%TAG2STR4   like int2str, but format integer (from 0 to 9999) with a constant length of 4
% 
% Prepared by Benoit Spinewine


if i<10
  string = ['000' int2str(i)];
elseif i<100
  string = ['00' int2str(i)];
elseif i<1000
  string = ['0' int2str(i)];
else
  string = int2str(i);
end;