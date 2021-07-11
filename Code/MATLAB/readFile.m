function [data,header]=readFile(file,headerBool,nocol)

 header=[];
 if nargin==1
   [nocol,header]=readHeader(file);
 end

 format='';
 for i=1:nocol
   format=strcat(format,'%f');
 end

 fid=fopen(file,'r');
 if nargin == 1 || headerBool
   textscan(fid,'%*s',nocol,'Delimiter',{' ','\t'},'CommentStyle','#','MultipleDelimsAsOne',1);
 end
 fileInCell=textscan(fid,format,'Delimiter',{'\b','\t',' '},'CommentStyle','%','MultipleDelimsAsOne',1);
 fclose(fid);

 data=zeros(length(fileInCell{1}),nocol);
 for i=1:nocol
   data(:,i)=fileInCell{i};
 end

end
