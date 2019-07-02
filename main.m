% Converts .ltcx to .inp with fully meshed hybrid beams and solid elements.
% Required: partsdefinition_workspace.mat
%           node.txt and connectivity.txt in ltcx format
% Define BCs and Mat-Props in CAE

%% Import lattice data as strut coordinates from modified .ltcx
% strutcoord x1 | y1 | z1 | x2 | y1 | z2 for a single strut on each row
clear;

[nodes, npath] = uigetfile('*.txt','Select NODE data to import');
[conn, cpath] = uigetfile('*.txt','Select CONNECTIVITY data to import');
N = [npath nodes];
C = [cpath conn];

% Extract nodes data: nodeid | x-coord | y-coord | z-coord
fid = fopen(N, 'rt');
node_data = textscan(fid, '%s id="%f" x="%f" y="%f" z="%f"/>');
fclose(fid);
nodes = [node_data{2} node_data{3} node_data{4} node_data{5}];

% Extract connectivity data: beamid | nodeid 1 | nodeid 2
fid2 = fopen(C, 'rt');
conn_data = textscan(fid2, '%s id="%f" n1="%f" n2="%f"/>');
fclose(fid2);
conn = [conn_data{2} conn_data{3} conn_data{4}];

% Rearrange to give matrix with: x1 | y1 | z1 | x2 | y1 | z2 on each row 
strutcoord = zeros(size(conn,1),6);

for i=1:size(conn,1)
    strutcoord(i,:)=[nodes(conn(i,2)+1, 2:4) nodes(conn(i,3)+1, 2:4)];
end
clear i nodes conn fid fid2 node_data conn_data N C; 

%% Subdivide strut into 3 sections and calculates length of each section
%        beam 20%           solid 60%            beam 20%
%     .------------.===========================.-------------.
% x1 y1 z1     x2 y2 z2                    x3 y3 z3       x4 y4 z4

strutlength = zeros(size(strutcoord,1),1);

% Determines the length of each strut
for i=1:size(strutcoord,1)
    strutlength(i,1) = sqrt(((strutcoord(i,4)-strutcoord(i,1))^2)+ ...
                   ((strutcoord(i,5)-strutcoord(i,2))^2)+ ...
                   ((strutcoord(i,6)-strutcoord(i,3))^2));
end

% Determine coordinates of each section
strutdir = zeros(size(strutcoord,1),3);
strut_sections = zeros(size(strutcoord,1),12);

for j=1:size(strutcoord,1)
    strutdir(j,:) = strutcoord(j,4:6)-strutcoord(j,1:3);
    strut_sections(j,1:3) = strutcoord(j,1:3);
    strut_sections(j,4:6) = (0.2.*strutdir(j,:))+strutcoord(j,1:3);
    strut_sections(j,7:9) = (0.8.*strutdir(j,:))+strutcoord(j,1:3);
    strut_sections(j,10:12) = strutcoord(j,4:6);
end

% Calculates length of each section and store in nX3 matrix
sectionlength = zeros(size(strutcoord,1),3);

for k=1:size(strutcoord,1)
    sectionlength(k,:) = [0.2*strutlength(k),0.6*strutlength(k),...
        0.2*strutlength(k)]; 
end

clear i j k;
%% Obtain translation and rotation of each section for assembly
% transformation: no. sections (struts x 3) X transformation
% x | y | z | xa | ya | za | xb | yb | zc | theta

% x, y, z corresponds to translations

% xa ya za xb yb zb theta corresponds to rotation as follows
% X-coordinate of point a on the axis of rotation.  xa
% Y-coordinate of point a on the axis of rotation.  ya
% Z-coordinate of point a on the axis of rotation.  za
% X-coordinate of point b on the axis of rotation.  xb
% Y-coordinate of point b on the axis of rotation.  yb
% Z-coordinate of point b on the axis of rotation.  zb
% Angle of rotation about the axis a–b, in degrees  theta

% Translation is applied before rotation

% a: [0 0 1], b: direction of section, c: axis of rotation  
% c = a x b = |a|.|b|.sin(theta).n
% theta = arccos((a.b)/(|a|.|b|))

% for i=1:size(strutcoord,1)
%     for j=1:3
%         transf((3*i)+(j-3),1:6) = [strut_sections(i,((3*j)-2):(3*j)),...
%             strut_sections(i,((3*j)-2):(3*j))];
%         transf((3*i)+(j-3),7:9) = cross(A,(strut_sections(i,((3*j)+1):...
%             ((3*j)+3))-strut_sections(i,((3*j)-2):(3*j))))+strut_sections...
%             (i,((3*j)-2):(3*j));
%     end
%     transf((3*i)-2:(3*i),10) = rad2deg(acos(dot(A,strutdir(i,:))/...
%         (norm(A)*norm(strutdir(i,:)))));
% end

transf = zeros(size(strutcoord,1)*3,10);
A = [0 0 10];

for i=1:size(strutcoord,1)
    for j=1:3
        transf((3*i)+(j-3),1:6) = [strut_sections(i,((3*j)-2):(3*j)),...
            strut_sections(i,((3*j)-2):(3*j))];
        transf((3*i)+(j-3),7:9) = cross(A,strutdir(i,:))+strut_sections...
            (i,((3*j)-2):(3*j));
    end
    transf((3*i)-2:(3*i),10) = rad2deg(acos(dot(A,strutdir(i,:))/...
        (norm(A)*norm(strutdir(i,:)))));
end

% ind = zeros(size(transf,1),3);
% for dfa=1:size(transf,1)
%     for dfr=1:3
%         ind(dfa,dfr) = ismember(transf(dfa,dfr+3),transf(dfa,dfr+6));
%     end
% end
% 
% para = zeros(size(coin,1),1);
% 
% for ds=1:size(ind,1)
%     para(ds)=~ismember(0,ind(ds,:));
% end
% 
% ind2 = find(para, sum(para(:) == 1));
% angles2=zeros(size(ind2,1),1);
% for oj=1:size(ind2,1)
%     angles2(oj) = transf(ind2(oj),10);
% end
% 
% angle2 = angle2';
% 



B = [10,0,0];

for k=0:size(strutdir,1)-1
    for m=1:3
        if transf((3*k)+m,10) == 180
            transf((3*k)+m,7:9) = cross(B,strutdir(k+1,:))+transf((3*k)+m,1:3);
        else
        end
    end
end


clear i j A;
%% Create parts with nodes and elements for all strut lengths
% snodes and bnodes: node no. | x | y | z
% selements and belements: element no. | 1st node | 2nd node | ...
% different strut lengths are stored in the third dimension
% usections: lengths of beam | solid | beam X unique sections

% Number of parts depend on number of unique strut lengths x 2
[usections, ~, ic] = unique(sectionlength,'rows','stable');
ustruts = size(usections,1);

load('partsdefinition_workspace3');

% Prompts user to select analysis type
anatype = questdlg('Please choose analysis type', ...
	'Analysis type', ...
	'Implicit','Explicit','Implicit');
switch anatype
    case 'Implicit'
        esolid = esolid_imp;
        ssolid = ssolid_imp;
        analysis = 'imp';
    case 'Explicit'
        esolid = esolid_exp;
        ssolid = ssolid_exp;
        analysis = 'exp';
end


% Create solid sections
snodes(:,:,ustruts) = zeros(size(ssolid,1),4);
selements(:,:,ustruts) = zeros(size(esolid,1),size(esolid,2));
for i=1:ustruts
   for j=1:size(ssolid,1)
       snodes(j,:,i) = [j, ssolid(j,1:2), ssolid(j,3).*usections(i,2)];
   end
   selements (:,:,i) = esolid;
end

% Create beam sections
bnodes(:,:,ustruts) = zeros(size(sbeam,1),4);
belements(:,:,ustruts) = zeros(1,4);
for m=1:ustruts
    for n=1:size(sbeam,1)
        bnodes(n,:,m) = [n, sbeam(n,1:3).*usections(m,1)];
    end
    belements (:,:,m) = ebeam;
end

clear ustruts i j m n ebeam esolid sbeam ssolid;

%% Create file

% Specify filename 
prompt = {'Filename','Job name','Model name'};
title = 'Input job details';
name = inputdlg(prompt,title,[1 50],{'Input file','Job-1','Model-1'});

% Selecting save directory
selpath = uigetdir('','Select save location');

% Create file. 'w' for overwrite permission, 'a' for append permission 
fid3=fopen([selpath '/' char(name{1}) '.inp'],'w');

%% Headings and Preprint

Headings = {'*Heading';
sprintf('** Job name: %s Model name: %s',name{2},name{3});
'** Generated by: Abaqus/CAE 6.14-4';
'*Preprint, echo=NO, model=NO, history=NO, contact=NO';
'**';
'** PARTS';
'**'};

for i=1:length(Headings)
    fprintf(fid3,'%s\n', Headings{i});
end
clear i;
%% Part. Loop for number of unique section lengths
% Solid parts
spn = permute(snodes,[2,1,3]);
spe = permute(selements, [2,1,3]);
for i=1:size(usections,1)
    
    % Part name
    fprintf(fid3,'*Part, name=Solid-Strut%u\n*Node\n',i);
    % Node
    fprintf(fid3,'      %u, %.9f, %.9f, %.9f\n',spn(:,:,i));
    switch anatype
        case 'Implicit'
            % Element
            fprintf(fid3,'*Element, type=C3D20R\n');
            fprintf(fid3,['%u, %u, %u, %u, %u, %u, %u, %u, %u, %u,'...
                '%u, %u, %u, %u, %u, %u,\n   %u, %u, %u, %u, %u\n'],...
                spe(:,:,i));
        case 'Explicit'
            % Element
            fprintf(fid3,'*Element, type=C3D8R\n');
            fprintf(fid3,['%u, %u, %u, %u, %u, %u, %u, %u, %u'...
                '\n'],spe(:,:,i));
    end
    % Nset
    fprintf(fid3,['*Nset, nset=Solid-Section-Nset%u, generate\n'...
        '1, %u, 1\n'],i, size(snodes,1));
    % Elset
    fprintf(fid3,['*Elset, elset=Solid-Section-Elset%u, generate\n'...
        '1, %u, 1\n'],i, size(selements,1));
    % Section
    fprintf(fid3,['** Section: Solid-Section\n*Solid Section, '...
        'elset=Solid-Section-Elset%u, material=Mat-Props-Solid\n,\n'],i);
    % End Part
    fprintf(fid3,'*End Part\n**\n');
end
% Wire parts
wpn = permute(bnodes,[2,1,3]);
wpe = permute(belements, [2,1,3]);
for j=1:size(usections,1)
    % Part name
    fprintf(fid3,'*Part, name=Wire-Strut%u\n*Node\n',j);
    % Node
    fprintf(fid3,'      %u, %.1f, %.1f, %.1f\n',wpn(:,:,j));
    % Element
    fprintf(fid3,'*Element, type=B32\n');
    fprintf(fid3, '%u, %u, %u, %u\n',wpe(:,:,j));
    % Nset
    fprintf(fid3,'*Nset, nset=Wire-Section-Nset%u, generate\n   1, 3, 1\n',j);
    % Elset
    fprintf(fid3,'*Elset, elset=Wire-Section-Elset%u\n   1\n',j);
    % Section
    fprintf(fid3,'** Section: Beam-Section Profile: Circ-Profile\n');
    fprintf(fid3,['*Beam Section, elset=Wire-Section-Elset%u, material=Mat-'...
        'Props-Beam, temperature=GRADIENTS, section=CIRC\n'],j);
    fprintf(fid3,'0.5\n0.124,-0.426,0.784\n');
    % End Part
    fprintf(fid3,'*End Part\n**\n');
end

clear i j wpe wpn;
%% Assembly
% Assembly name
fprintf(fid3,'**\n** ASSEMBLY\n**\n*Assembly, name=Assembly\n**\n');

% Solid part
for i=2:3:(size(transf,1)-1)
    fprintf(fid3,'*Instance, name=Solid-Strut-%u, part=Solid-Strut%u\n',...
        (i+1)/3, ic((i+1)/3));
    fprintf(fid3,'          %.9f, %.9f, %.9f\n', transf(i,1:3));
    if transf(i,10) == 0
    else
        fprintf(fid3,['          %.9f, %.9f, %.9f, %.9f, %.9f, %.9f,'...
            '%.3f\n'], transf(i,4:10));
    end
    fprintf(fid3,'*End Instance\n**\n');
end

% Wire part
index = zeros(size(transf,1)*2/3,1);
for j=1:size(index,1)/2
    index((2*j)-1) = (3*j)-2;
    index(2*j) = (3*j);
end

labo=1;
labe=1;
for k=1:size(index,1)
    if mod(k,2) == 0
        fprintf(fid3,'*Instance, name=Wire-Strut-%u%u, part=Wire-Strut%u\n',...
            labe, 2, ic(index(k)/3));
        labe=labe+1;
    else
        fprintf(fid3,'*Instance, name=Wire-Strut-%u%u, part=Wire-Strut%u\n',...
            labo, 1, ic((index(k)+2)/3));
        labo=labo+1;
    end
    fprintf(fid3,'          %.9f, %.9f, %.9f\n', transf(index(k),1:3));
    if transf(index(k),10) == 0
    else
        fprintf(fid3,['          %.9f, %.9f, %.9f, %.9f, %.9f, %.9f,'...
            '%.3f\n'], transf(index(k),4:10));
    end
    fprintf(fid3,'*End Instance\n**\n');
end

% Create Nsets for constraints
% Both ends of beam to be assigned to an Nset
% One end for coupling constraints at beam-solid interface
% Other end for tie constraints at beam junctions
% All slave nodes at the same junction in the same Nset

% Value in odd rows of 'index' corresponds to beams in 'transf' where end 1
% is at the junction and end 2 is at the interface (see 'strutcoord')
% Value in even rows of 'index' corresponds to beams in 'transf' where end
% 3 is at the interface and end 4 is at the junction (see 'strutcoord')

% Create Nset for nodes on beam at interface
ei = 0;
oi = 0;
for ni=1:size(index,1)
    if mod(ni,2) == 0               % even rows
        ei=ei+1;
        fprintf(fid3,['*Nset, nset=WS-%u%u-NsetInt-%u, instance=Wire-Strut'...
            '-%u%u\n %u,\n'], ei, 2, ni, ei, 2, 1);
    else                            % odd rows
        oi=oi+1;
        fprintf(fid3,['*Nset, nset=WS-%u%u-NsetInt-%u, instance=Wire-Strut'...
            '-%u%u\n %u,\n'], oi, 1, ni, oi, 1, 3);
    end
end

% Create Nset for nodes at junctions

% Create an juncnodes (junction nodes) of   x1 y1 z1
%                                           x4 y4 z4    from strut_sections

juncnodes = zeros(size(index,1),3);
for in=1:2:size(juncnodes,1)-1
    juncnodes(in,:) = strut_sections((in+1)/2,1:3);
    juncnodes(in+1,:) = strut_sections((in+1)/2,10:12);
end

% Loop through all nodes and indentify nodes that are near each other in 
% juncnodes and obtain indexS
% if index of loop is in indexS skip index, if not then store in indexM
% find all cells where rs = value and remove

indexM = zeros(size(juncnodes,1),1);
indexRS = cell(size(juncnodes,1),1);
indexS = zeros(size(juncnodes,1),20);        % Up to 20 slave nodes

for rs=1:size(juncnodes,1)
    if ismember(rs,indexS) == 1
    else
        indexM(rs)= rs;
        indexRS(rs) = rangesearch(juncnodes(:,:),juncnodes(rs,:),0.01);
        indexS(rs, 1:cellfun('length',indexRS(rs))) = cell2mat(indexRS(rs));
        indexS(rs, find(indexS(rs,:)==rs))=0; %#ok<FNDSB> % remove Mnode
    end
end

% Create Nset of masters (intnodes(indexM)) and Nset of slaves
% if no slave exists, do not create master
% if index is even, instance = Wire-Strut-('index'/2)2
% if index is odd, instance = Wire-Strut-(('index'+1)/2)1
s = 0;    %count number of surfaces to create

for im=1:size(indexM,1)
    if indexM(im) ~= 0 && ~isempty(find(indexS(im,:))) %#ok<EFIND>        
        if mod(indexM(im),2) == 0   % Master, even rows
            s = s+1;
            fprintf(fid3,['*Nset, nset=WS-%u%u-NsetJunM-%u, instance='...
                'Wire-Strut-%u%u\n %u,\n'], indexM(im)/2, 2, s,...
                indexM(im)/2, 2, 3);
            ms(s,:) = [indexM(im)/2, 2, ];
            for is=1:size(indexS,2) % Slaves
                if indexS(im,is) == 0
                elseif indexS(im,is) ~= 0 && mod(indexS(im,is),2) == 0
                    fprintf(fid3,['*Nset, nset=WS-%u%u-NsetJunS-%u, '...
                        'instance=Wire-Strut-%u%u\n %u,\n'],...
                        indexM(im)/2, 2, s, indexS(im,is)/2, 2, 3);
                else
                    fprintf(fid3,['*Nset, nset=WS-%u%u-NsetJunS-%u, '...
                        'instance=Wire-Strut-%u%u\n %u,\n'],...
                        indexM(im)/2, 2, s, (indexS(im,is)+1)/2, 1, 1);
                end
                clear is
            end
        else                        % Master, odd rows
            s = s+1;
            fprintf(fid3,['*Nset, nset=WS-%u%u-NsetJunM-%u, instance='...
                'Wire-Strut-%u%u\n %u,\n'], (indexM(im)+1)/2, 1, s,...
                (indexM(im)+1)/2, 1, 1);
            ms(s,:) = [(indexM(im)+1)/2, 1]; %#ok<*SAGROW>
            for is=1:size(indexS,2) % Slaves
               if indexS(im,is) == 0
               elseif indexS(im,is) ~= 0 && mod(indexS(im,is),2) == 0
                   fprintf(fid3,['*Nset, nset=WS-%u%u-NsetJunS-%u, '...
                       'instance=Wire-Strut-%u%u\n %u,\n'],...
                       (indexM(im)+1)/2, 1, s, indexS(im,is)/2, 2, 3);
               else
                   fprintf(fid3,['*Nset, nset=WS-%u%u-NsetJunS-%u, '...
                       'instance=Wire-Strut-%u%u\n %u,\n'],...
                       (indexM(im)+1)/2, 1, s, (indexS(im,is)+1)/2, 1, 1);
               end
               clear is
            end
        end
    else
    end
end

% Create Elset for both ends of each solid instance
% For each Elset, create a Surface for coupling constraint
switch anatype
    case 'Implicit'
        elset = {'4,  20,  4';'1,  17,  4'};
        for o=2:3:(size(transf,1)-1)    % loop through each solid instance
            for p=2:-1:1                    % two surfaces
                fprintf(fid3,['*Elset, elset=SS-%u-Elset-S%u, instance=Solid-Strut'...
                    '-%u, generate\n  ' elset{p} '\n'], (o+1)/3, p, (o+1)/3);
                fprintf(fid3,['*Surface, type=ELEMENT, name=SS-%u-Surf-S%u\n'...
                    'SS-%u-Elset-S%u, S%u\n'], (o+1)/3, p, (o+1)/3, p, p);
            end
        end
    case 'Explicit'
        elset = {'1,  2,  3, 13, 14, 15, 25, 26, 27, 37, 38, 39';...
            '10, 11, 12, 22, 23, 24, 34, 35, 36, 46, 47, 48'};
        surf = {1; 2};
        for o=2:3:(size(transf,1)-1)    % loop through each solid instance
            for p=2:-1:1                    % two surfaces
                fprintf(fid3,['*Elset, elset=SS-%u-Elset-S%u, instance=Solid-Strut'...
                    '-%u\n  ' elset{p} '\n'], (o+1)/3, surf{p}, (o+1)/3);
                fprintf(fid3,['*Surface, type=ELEMENT, name=SS-%u-Surf-S%u\n'...
                    'SS-%u-Elset-S%u, S%u\n'], (o+1)/3, surf{p}, (o+1)/3,...
                    surf{p}, surf{p});
            end
        end
end


% Create surfaces for beam junction
% One for master and other for slave

% for s=1:s(end)
%     fprintf(fid3,['*Surface, type=Node, name=WS-%u-Surf-M\n'...
%         'WS-%u%u-NSetJunM-%u, 1\n'], s, ms(s,1), ms(s,2), s);
%     fprintf(fid3,['*Surface, type=Node, name=WS-%u-Surf-S\n'...
%         'WS-%u%u-NSetJunS-%u, 1\n'], s, ms(s,1), ms(s,2), s);
%     
% end

% Create tie constraint (beam junction)

% for s=1:s(end)
%     fprintf(fid3,['** Constraint: N-N-Tie-%u\n*Tie, name=N-N-Tie-%u, '...
%         'adjust=no\nWS-%u-Surf-S, WS-%u-Surf-M\n'], s, s, s, s);
% end

% Create mpc constraint (beam junction)

for s=1:s(end)
    fprintf(fid3,['** Constraint: N-N-mpc-%u\n*MPC\nTIE, '...
        'WS-%u%u-NSetJunS-%u, WS-%u%u-NSetJunM-%u\n'], s, ms(s,1),...
        ms(s,2), s, ms(s,1), ms(s,2), s);
end


% Create coupling constraint (solid to beam)

rn = ei+oi;

switch anatype
    case 'Implicit'
        for n=1:rn(end)
            if mod(n,2) == 0
                fprintf(fid3,['** Constraint: N-S-%u\n*Coupling, constraint name='...
                    'N-S-%u, ref node=WS-%u%u-NsetInt-%u, surface=SS-%u-Surf-S%u\n'...
                    '*Distributing, weighting method=UNIFORM\n'], n, n, n/2, ...
                    (mod(n,2)*-1)+2, n, n/2, mod(n,2)+1);
            else
                fprintf(fid3,['** Constraint: N-S-%u\n*Coupling, constraint name='...
                    'N-S-%u, ref node=WS-%u%u-NsetInt-%u, surface=SS-%u-Surf-S%u\n'...
                    '*Distributing, weighting method=UNIFORM\n'], n, n, (n/2)+0.5, ...
                    (mod(n,2)*(-1))+2, n, (n/2)+0.5, mod(n,2)+1);
            end
        end
    case 'Explicit'
        for n=1:rn(end)
            if mod(n,2) == 0
                fprintf(fid3,['** Constraint: N-S-%u\n*Coupling, constraint name='...
                    'N-S-%u, ref node=WS-%u%u-NsetInt-%u, surface=SS-%u-Surf-S%u\n'...
                    '*Distributing, weighting method=UNIFORM\n'], n, n, n/2, ...
                    (mod(n,2)*-1)+2, n, n/2, mod(n,2)+2);
            else
                fprintf(fid3,['** Constraint: N-S-%u\n*Coupling, constraint name='...
                    'N-S-%u, ref node=WS-%u%u-NsetInt-%u, surface=SS-%u-Surf-S%u\n'...
                    '*Distributing, weighting method=UNIFORM\n'], n, n, (n/2)+0.5, ...
                    (mod(n,2)*(-1))+2, n, (n/2)+0.5, mod(n,2));
            end
        end        
end

fprintf(fid3,'*End Assembly\n**\n');

clear labo labe i j k m n o p elset ei oi ni is im in rs s rn
%% Materials

Materials = {'** MATERIALS';
    '**';
    '*Material, name=Mat-Props-Solid';
    '*Density';
    ' 8e-06,';
    '*Elastic';
    '1930., 0.3';
    '*Plastic';
    '170., 0.';
    '175.,10.'};

for i=1:length(Materials)
    fprintf(fid3,'%s\n', Materials{i});
end

Materials = {'*Material, name=Mat-Props-Beam';
    '*Density';
    ' 8e-06,';
    '*Elastic';
    '1930., 0.3';
    '*Plastic';
    '170., 0.';
    '175.,10.';
    '**'};

for j=1:length(Materials)
    fprintf(fid3,'%s\n', Materials{j});
end
clear i j;

%% Step

switch anatype
    case 'Implicit'
        load = {'** ----------------------------------------------------------------';
            '** ';
            '** STEP: Load';
            '**';
            '*Step, name=Load, nlgeom=YES, inc=1000';
            '*Static';
            '0.1, 1., 1e-05, 1.';
            '**';
            '** OUTPUT REQUESTS';
            '** ';
            '*Restart, write, frequency=0';
            '** ';
            '** FIELD OUTPUT: F-Output-1';
            '**';
            '*Output, field';
            '*Node Output';
            'CF, RF, U';
            '*Element Output, directions=YES';
            'LE, PE, PEEQ, PEMAG, S, SF';
            '**';
            '** HISTORY OUTPUT: H-Output-1';
            '**';
            '*Output, history, variable=PRESELECT';
            '*End Step'};
    case 'Explicit'
        load = {'** ----------------------------------------------------------------';
            '** ';
            '** STEP: Load';
            '**';
            '*Step, name=load, nlgeom=YES';
            '*Dynamic, Explicit';
            ', 5.';
            '*Bulk Viscosity';
            '0.06, 1.2';
            '**';
            '** OUTPUT REQUESTS';
            '** ';
            '*Restart, write, number interval=1, time marks=NO';
            '** ';
            '** FIELD OUTPUT: F-Output-1';
            '**';
            '*Output, field';
            '*Node Output';
            'A, RF, U, V';
            '*Element Output, directions=YES';
            'LE, PE, PEEQ, S, SF';
            '**';
            '** HISTORY OUTPUT: H-Output-1';
            '**';
            '*Output, history, variable=PRESELECT';
            '*End Step'};   
end

for i=1:length(load)
    fprintf(fid3,'%s\n', load{i});
end
clear i;

%% Complete
CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'modal';
complete = msgbox({'\fontsize{14} File created successfully';...
    ' \fontsize{12}';' Import into Abaqus/CAE to define:';...
    '     1. Material properties';'     2. Steps';...
    '     3. Boundary Conditions';'     4. Output requests'},...
    'Complete',CreateStruct);
