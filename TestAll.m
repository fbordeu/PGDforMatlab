function res = TestAll(d)
%Function to test the integrity of the library
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%

%% this make the script runnable for any directory
cpw = pwd();
libpath = mfilename('fullpath');
cd(libpath(1:find(libpath==filesep(),1,'last')))

res = true;


logfid = fopen('TestAll.log','w');

fprintf(logfid,'Start Date %s \n',date);
c = clock;
fprintf(logfid,'Start Time %ih%i:%s \n',c(4:5),num2str(floor(c(6))));

if ~exist('d','var')
 content = dir();
 disp('0) Execute all tests')
 cpt = 1;
 testnames ={};
 for i = 1:numel(content)
     d = content(i);
     if d.name(1) == '.'
         continue
     end
     if d.isdir
        a = dir([d.name filesep() 'Tests' filesep() '*test*m']);
        a = [a ; dir([d.name filesep() 'Tests' filesep() '*test*p']) ];
        if numel(a)>0
         disp([num2str(cpt) ') ' d.name])
         testnames{cpt} = d.name;
         cpt = cpt +1;
        end           
     end
 end

 s=input('Please choose the toolbox to test (enter to test all): ' ,'s');
 switch s
     case {'0',''}
         d = 'All';
     otherwise
         d = testnames{str2double(s)}; 
 end
end
        
p = pwd;
if exist('d','var')
    if strcmpi('All',d)
        content = dir();    
    else
       content = dir(['*' d]);
    end
else
    content = dir();
end

for i = 1:numel(content)
    d = content(i);
    if d.name(1) == '.'
        continue
    end
    if d.isdir
        fprintf(1,' Working in directory %s **********************************\n', d.name);
        fprintf(logfid,' Working in directory %s *************************************\n', d.name);
        cd(p);
        res = res & TestDirectory(p,d.name,logfid);
    end
    
end
c = clock;
fprintf(logfid,'Stop Date %s %ih%i:%s \n',date,c(4:5),num2str(floor(c(6))));
fclose(logfid);
cd(p);

cd(cpw);
end


function status = TestDirectory(p,d,logfid)

status = true;
cd([p filesep d ]) 

a = dir(['Tests' filesep '*test*m']);
a = [a ; dir(['Tests' filesep '*test*p'])];


fprintf(logfid, ' Found %i test\n', numel(a) );
cp = pwd();
outputstring =  {'NOT OK '  'OK'};
for i = 1:numel(a)
    f = a(i);
    cd(cp);
    cd('Tests')
    fprintf(logfid,'  Function %s \n',f.name);
    
    c = clock;
    fprintf(logfid,'   Start Date %s %ih%i:%s \n',date,c(4:5),num2str(floor(c(6))));
    try
        disp(f.name(1:end-2))
        istatus =  feval(f.name(1:end-2));
    catch ME
        ME.identifier;
        fprintf(logfid,'   ERROR :%s \n',ME.message);
        fprintf(logfid,'   Stack :\n');
        for i =1:numel(ME.stack)
            fprintf(logfid,'    %s\n', ['Function ' ME.stack(i).name ' ' ME.stack(i).file ':' num2str(ME.stack(i).line)]);
        end
        istatus = false;
    end
    c = clock;
    fprintf(logfid,'   Stop Date %s %ih%i:%s \n',date,c(4:5),num2str(floor(c(6))));
    fprintf(logfid,'   Status %s \n',outputstring{istatus+1});
    fprintf(2-istatus,'   Status %s \n',outputstring{istatus+1});
    status = status && istatus;
end



end