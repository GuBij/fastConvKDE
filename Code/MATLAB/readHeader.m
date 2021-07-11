function [nocol,header]=readHeader(file)

 fid=fopen(file,'r');
 header=textscan(fid,'%[abcdefghijklmnopqrstuvwxyzACKUVPIDMRS0123456789_:-/]','Delimiter',{',','"','\t','\b'},'CommentStyle','#','MultipleDelimsAsOne',1);
 fclose(fid);
 header=header{1};

 nocol=length(header);

 count=0;
 while ~isempty(str2num(header{nocol-count}))
  count=count+1;
 end
 nocol=nocol-count;

end

