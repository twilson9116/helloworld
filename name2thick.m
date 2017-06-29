for i=2:length(names)
str=names{i};
cut=strfind(str,'-');
cut=cut(1)-1;
str=str(1:cut);
str=double(str);
str=str(str>47 & str<58);
thickness(i)=str2num(char(str));

end


thickness=thickness';