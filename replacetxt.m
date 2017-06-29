function y=replacetxt(file);
fin = fopen(file,'w');
asdfasdf
while ~feof(fin)
   s = fgetl(fin)
   s = strrep(s,'.01000','.00000');
   s = strrep(s,'.02000','.00000');
   s = strrep(s,'.03000','.00000');
   fprintf(fout,'%s\n',s);
   %disp(s)
end
fclose(fin)

end