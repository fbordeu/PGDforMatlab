function res = Kompressor1test()
%
% This file is subject to the terms and conditions defined in
% file 'LICENSE.txt', which is part of this source code package.
%
% Principal developer : Felipe Bordeu (Felipe.Bordeu@ec-nantes.fr)
%

disp('Running Kompressor1test')

plott = true;
res = true;

x = (0:2800)/2800;
y = (0:2700)/2700;
z = (0:2600)/2600;
a = (0:2500)/2500;
b = (0:2400)/2400;


FF = cell(2,1);
FF{1} = [x'*0.1 (x.^2)'*100 x'*1000 x'*10000    sin(x*3)'*10000 cos(x)'*10 ];
FF{2} = [y' y'   (y.^3)'  y' sin(y)' cos(y)'   ];
FF{3} = [z' z' z' (z.^4)'    sin(z)' cos(z)'  ];
%FF{4} = [(a.^2)' a' a' (a.^4)'    sin(a)' cos(a)'  ];
%FF{5} = [b' (b.^2)' b' (b.^2)'    sin(10*b)' sin(b)'  ];

for i = 1:numel(FF)
    FF{i} = repmat(FF{i},1, 4);
end

tic
disp('Solving with mex')
sol = recompact(FF,'verbose',true,'max_added_modes',6);
ctime = toc;
tic
disp('Solving with V6')
bb = recompact(FF,'usev6',true,'verbose',true,'max_added_modes',6);
mtime = toc;

disp(['ctime  ' num2str(ctime)])
disp(['mtime  ' num2str(mtime)])

NumberOfdims = size(FF,1);

if plott; 
    figure ; %#ok<*UNRCH>
end 
disp(['c modes '   num2str(size(sol{1},2))])
disp(['M modes '   num2str(size(bb{1},2))])


nmodes = min(size(sol{1},2), size(bb{1},2));
cpt =1; %#ok<*NASGU>
for mode = 1:nmodes
 for i =1:NumberOfdims 
    d = dot(bb{i}(:,mode),sol{i}(:,mode))/(norm(bb{i}(:,mode))*norm(sol{i}(:,mode)) )  ;
    
    if plott; 
      subplot(nmodes,NumberOfdims,cpt)
      cpt = cpt +1;
      plot(bb{i}(:,mode),'m'); 
      hold on;
      plot(sol{i}(:,mode)*sign(d),'c--') 
      title(num2str(d  ));
    end
    
    if abs(d) < 0.999
       res = false;
    end
    

  end
end


%pnet(con,'close');
disp('DONE Kompressor1test')
return 




dataint = readpxdmf('/home/fbordeu/Dropbox/ParaView Data2/ParaMat3D.pxdmf')
stime = tic; data = recompact(dataint,'useC',true); ctime = toc(stime)
stime = tic; data = recompact(dataint); mtime = toc(stime)

load FF;
NumberOfdims = size(FF,1);
pnet_putvar(con,NumberOfdims);
for i=1:size(FF,1)
    %pause(1)
    b = FF{i,1}';
    pnet_putvar(con,b);
end
%pnet(con,'close');
%return 

sol= {};
for i=1:size(FF,1)
    pause(1)
    sol{i} =  pnet_getvar(con);
end
% 
ff{1} = FF{1}';
ff{2} = FF{2}';
ff{3} = FF{3}';
ff{4} = FF{4}';
ff{5} = FF{5}';
tic; bb = recompact(ff'); toc

figure ;
for i =1:5
    subplot(2,5,i)
    plot(bb{i}(:,1),'m'); hold on;
    subplot(2,5,i+5)
    plot(sol{i}(:,1),'c')
end

figure ;
for i =1:5
    subplot(2,5,i)
    plot(bb{i}(:,2),'m'); hold on;
    subplot(2,5,i+5)
    plot(sol{i}(:,2),'c')
end

figure ;
for i =1:5
    subplot(2,5,i)
    plot(bb{i}(:,3),'m'); hold on;
    subplot(2,5,i+5)
    plot(sol{i}(:,3),'c')
end

figure ;
for i =1:5
    subplot(2,5,i)
    plot(bb{i}(:,6),'m'); hold on;
    subplot(2,5,i+5)
    plot(sol{i}(:,6),'c')
end


%figure ;
%plot(bb{2}(:,2),'m'); hold on;
%plot(sol{2}(:,2),'c')


pnet(con,'close');

%return 