fileID = fopen('TypListe.csv','r');
ListeArray = textscan(fileID,'%s%u%u%u%f');
fclose(fileID);
ListeArray{5}


fid = fopen('TypListe2.csv','w');fprintf(fid,'"%s" %i %i %i %.2f',col,rgb(1),rgb(2),rgb(3),rnd);fclose(fid);

%%
Name='Mond';
Farbe='Grau normal';
R1=85; G1=86; B1=87;
R = R1/255; G = G1/255; B = B1/255;
Hex='555657';
fid = fopen('PlanetenFarben2.csv','a');
fprintf(fid,'%.4f;%.4f;%.4f;%s;%s;%s;%i;%i;%i\r',R,G,B,Name,Farbe,Hex,R1,G1,B1);
fclose(fid);

