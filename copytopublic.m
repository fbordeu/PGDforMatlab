function copytopublic(publiclibname)
restoredefaultpath();
addpath(pwd);
setMatlabTools('nosave');
assert(strcmp([pwd filesep mfilename], mfilename('fullpath')),'Error copytopublic cannot be launched from another directory');
key=123456;
copytopublic_priv(publiclibname,key);
end

function copytopublic_priv(publiclibname,key)
%% Dont run this scrip inside matlab
%% this file is called by the script file distribute

%% we use the key to be sure that you dont use this file manualy 

try
    if key ~= 123456
        assert()
    end
catch
    publiclibname
    disp('ERROR!! : Run this script from the distribute.sh file. Don''t Forget to config the destination path')
    return 
end

% script to pcodedify all the matlab tools
clc

%recover the current path
cpw = pwd();


% recover the path of the library
libpath = mfilename('fullpath');
thisfile = [libpath(find(libpath==filesep(),1,'last')+1:end) '.m'];
libpath = libpath(1:find(libpath==filesep(),1,'last'));
% build the path of the public library
publiclibpath = [libpath(1:find(libpath==filesep(),1,'last')) '..' filesep() '..' filesep() publiclibname filesep() 'MatlabTools' filesep() ]

% 
%CopyOnlyDoc('/home/fbordeu/projects/PGD/MatlabTools/setMatlabTools.m', '/home/fbordeu/projects/PGDPublic/MatlabTools/setMatlabTools.m')
%return



%%
cd(libpath)

% if the lib directory is dirty no distibution is made
[~,out] = system('svn diff');
if numel(out) > 0
    error('The source directory is dirty no distibution possible');
    return
end

disp('Reatriving list form svn ...')

[~,filesInSvn] = system('svn list --recursive');

filesInSvn = cellfun(@(arg)( [ libpath arg] ),strsplit(filesInSvn),'UniformOutput', false);

exttocopy = {'.geo' '.msh' '.txt'};
filestoskip = { thisfile };

%copy the tree structure
GeneratePcode(libpath, publiclibpath,exttocopy,filestoskip,filesInSvn,0);

%% to generate a version file
[~,result] = system('svn info --xml '); 
data1= strfind(result,'vision'); 
data2 = strfind(result(data1:end),'"'); 
version = (result((data1+data2(1)):(data1(1)+data2(2)-2)));

fid = fopen( [publiclibpath 'version'],'w');
fprintf(fid,'%s\n',version);
fclose(fid);

%%
cd(cpw)

end

function GeneratePcode(libpath,publibpath,exttocopy,filestoskip,filesInSvn,indent)

% only copy file with this extentions

cpw = pwd();
cd(libpath)
content = dir() ;
for i = 1:numel(content)
    d = content(i);
    if d.name(1) == '.'
        continue
    end
    sep = repmat('  ',1,indent);
    
    filename = [libpath d.name];
    %&& exist([d.name filesep 'Src' filesep],'dir')
    if d.isdir 
        % mkdir for all the directories
        disp ([sep d.name ])
        [status,message,messageid] = mkdir(publibpath,d.name);
        %% if the path contain 'data' we copy also the pxdmf files
        if (strfind(d.name,'data'))
            GeneratePcode([filename filesep()] ,[publibpath d.name filesep() ],{exttocopy{:} '.pxdmf' '.bin'},filestoskip,filesInSvn,indent+1)
        else
            GeneratePcode([filename filesep()] ,[publibpath d.name filesep() ],exttocopy,filestoskip,filesInSvn,indent+1)
        end
        
    else
        [~,name,ext] = fileparts(d.name) ;
      if sum(cellfun(@(arg) strcmpi(arg,d.name), filestoskip)) == 0 
        if sum(cellfun(@(arg) strcmpi(arg,ext), exttocopy)) >= 1 
            if CopyFileFromSVN(filename,[publibpath d.name],filesInSvn )
                disp ([sep d.name ])    
            end
        else
            if strcmpi('.m',ext) 
                
                if numel(strfind(lower(filename), 'example')) > 0
                    if CopyFileFromSVN(filename,[publibpath d.name],filesInSvn )
                        disp ([sep d.name ])
                    end
                else
                    cd(publibpath)
                    if PcodeFileFromSVN(filename,filesInSvn)
                        disp ([sep '(Pcoded) ' d.name ' (doc)'])
                    end
                    cd(cpw)
                end
            else
                if strfind(ext,'.mex') 
                    if CopyFileFromSVN(filename,[publibpath d.name],filesInSvn )
                        disp ([sep d.name ])
                    end
                end
                
            end
            
        end
      else
         disp ([sep '(Skip)   ' d.name ])
      end

        
    end

end
cd(cpw);
end

function res =CopyFileFromSVN(from,to,filesInSvn )
  res = false;
  if sum(cellfun(@(arg) strcmpi(arg,from), filesInSvn)) >= 1 
      copyfile(from,to)
      res = true;
  end  
end

function res = PcodeFileFromSVN(filename,filesInSvn )
  res = false;
  if sum(cellfun(@(arg) strcmpi(arg,filename), filesInSvn)) >= 1 
      %setenv('INFILE',filename); 
      %setenv('OUTFILE',[filename(find(filename==filesep(),1,'last')+1:end) ]); 
      %!perl -wne '/^%/?$a=1&&print:$a?last:0' $INFILE > $OUTFILE

      % to delete the old compiled files
      last = strfind(filename,filesep());
      if(exist([filename((last(end)+1):end-1) 'p'],'file'))
          delete([filename((last(end)+1):end-1) 'p'])
      end
      GenerateDoc(filename,filename(find(filename==filesep(),1,'last')+1:end) )

      %CopyOnlyDoc(filename,filename(find(filename==filesep(),1,'last')+1:end) )
      pcode(filename)
      res = true;
  end  
end

function GenerateDoc(fin,fout)
 [~,~,~] = mkdir('doc');
 fout = [ 'doc' filesep() fout];
 FunctionOrClassName  = fin(find(fin==filesep(),1,'last')+1:end-2);
 htmlstr = help2html(FunctionOrClassName,'','-doc');
 
 if exist(FunctionOrClassName,'class') == 8
    htmlstr = changeLinks(htmlstr,FunctionOrClassName);
 end
 
 fid = fopen([fout(1:end-2) '.html'],'w');

 fprintf(fid,'%s',htmlstr);
 fclose(fid);
 if exist(FunctionOrClassName,'class') == 8
     [~,~,~] = mkdir(['doc' filesep()  FunctionOrClassName])
     out = methods(FunctionOrClassName);
     for i = 1:numel(out)
         htmlstr = help2html([FunctionOrClassName '.' out{i}],'','-doc');
         fid = fopen(['doc' filesep() FunctionOrClassName filesep() out{i} '.html'],'w');
         fprintf(fid,'%s',htmlstr);
         fclose(fid);
     end
 end
end

function output = changeLinks(input,classname)

% we start from the last 
patter1 = 'href="';
n = numel(patter1);
indexs = fliplr(strfind(input,patter1));

[M1,M2] = methods(classname);
for i = indexs
   % we recover the name of the function    
   i1 = strfind(input(i:end),'.');
   i2 = strfind(input(i:end),''')');
   functionname = input((i+i1):(i+i2-2));
   res = cellfun(@(arg1)(strcmp(arg1, functionname)),M2);
   resstr = [M2{res(:,3),4} 'matlab'];
   %keyboard
   if strcmp(resstr(1:3),'mat') || strcmp(resstr(1:3),'han')
       continue
   end
   input = strrep(input, input((i+6):i+i2), [classname '/' functionname '.html']);
end
output = input;

%matlab:doc('ELEMENT.hammer_points')
%"ELEMENT/hammer_points"


end

function CopyOnlyDoc(fin, fout) %#ok<DEFNU>
  fiid = fopen(fin,'r');
  foid = fopen(fout,'w');
  tline = fgets(fiid);
  incopy = numel(strfind(lower(tline),'%'))>0 ;
  while ischar(tline)
      
    S = strtrim(tline);
    if incopy 
        if numel(S)>0 && strcmp(S(1),'%')
            fprintf(foid,'%s ',tline);
        else 
            incopy = false;
        end
    else
        if numel(strfind(lower(S),'function'))>0 
            incopy = true;
            fprintf(foid,'%s ',tline);
        end
    end
    tline = fgets(fiid);
  end
  fclose(fiid);
  fclose(foid);
end
