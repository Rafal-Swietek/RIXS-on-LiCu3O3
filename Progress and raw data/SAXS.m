function SAXS
%
% MATLAB function from:
% 'Elements of Modern X-ray Physics' by Jens Als-Nielsen and Des McMorrow
%
% Calculates: Calculates SAXS from a sphere, disk and rod
% Calls to:
close all; clear all;

Q=0:0.001:2; Q2=0.1:0.002:1; Q3=0.11:0.002:1; Q4=0.25:0.002:2;
FS=12; % Font size
disp('Running: may take some time to complete')

subplot(1,2,1) % Plot on linear scale
axlm=[0 10 0 1.1]; axis(axlm);
c1=[0.9 0.9 1];
patch([axlm(1) axlm(2) axlm(2) axlm(1)],[axlm(3) axlm(3) axlm(4) axlm(4)],[0 0 0 0],c1);
set(gca,'layer','top')

si = @(x) sinint(x); %sine integral function
% Sphere

R=50;
Rg=sqrt(3/5)*R;
F1=3*(sin(Q*R)-(Q*R).*cos(Q*R))./(Q*R).^3;
hs=line(Q*Rg,F1.*F1,'color','b','linewidth',1.5,'linestyle',':');

% Disk

R=50;
Rg=sqrt(1/2)*R;
p17=2./(Q*R).^2.*(1-besselj(1,2*Q*R)./(Q*R));
hd=line(Q*Rg,p17,'color','r','linewidth',1.5,'linestyle','--');

% Rod

L=50;
for ii=1:length(Q)
    x=Q(ii)*L;
    p15(ii)=2*si(x)/x-4*sin(x/2).^2./(x)^2;
end
%quadl(si,0,x)
Rg=sqrt(1/12)*L;
hr=line(Q*Rg,p15,'color','m','linewidth',1.5,'linestyle','-.');

set(gca,'xscale','linear','yscale','linear','fontname','times',...
'fontsize',FS,'linewidth',1.0,'gridlinestyle',':')
ylabel('$$\left| \mathcal F(\mathrm Q)\right|\,^2$$','interpreter','latex')
text(0.5,-0.1,'$$\mathrm Q R_g$$','interpreter','latex','FontName','Times',...
'FontSize',12,'horizontalalignment','center','units','normalized')
axis square
box on; grid on
legend([hs hd hr],'Sphere','Disk','Rod')
text(0.07,0.07,'(a)','FontName','Times','FontSize',12,...
'interpreter','latex','units','normalized')

subplot(1,2,2) % Plot on Log scale
axlm=[1 20 1e-4 1.1]; axis(axlm);
c1=[0.9 0.9 1];
patch([axlm(1) axlm(2) axlm(2) axlm(1)],[axlm(3) axlm(3) axlm(4) axlm(4)],...
[0 0 0 0],c1);
set(gca,'layer','top')

% Sphere
R=50;
Rg=sqrt(3/5)*R;
F1=3*(sin(Q*R)-(Q*R).*cos(Q*R))./(Q*R).^3;
hs=line(Q*Rg,F1.*F1,'color','b','linewidth',1.5,'linestyle',':');
line(Q2*Rg,0.0000017./Q2.^4,'color','b','linewidth',1.0,'linestyle','-')

% Disk

R=50;
Rg=sqrt(1/2)*R;
p17=2./(Q*R).^2.*(1-besselj(1,2*Q*R)./(Q*R));
hd=line(Q*Rg,p17,'color','r','linewidth',1.5,'linestyle','--');
line(Q3*Rg,0.00095./Q3.^2,'color','r','linewidth',1.0,'linestyle','-')

% Rod

L=50;
for ii=1:length(Q)
    x=Q(ii)*L;
    p15(ii)=2*si(x)/x-4*sin(x/2).^2./(x)^2;
end

Rg=sqrt(1/12)*L;
hr=line(Q*Rg,p15,'color','m','linewidth',1.5,'linestyle','-.');
line(Q4*Rg,0.075./Q4.^1,'color','m','linewidth',1.0,'linestyle','-')

set(gca,'xscale','log','yscale','log','fontname','times','fontsize',FS,...
'linewidth',1.0,'minorgridlinestyle','none','gridlinestyle',':')
text(0.07,0.07,'(b)','FontName','Times','FontSize',12,...
'interpreter','latex','units','normalized')
text(0.7,0.28,'$$1/\mathrm Q^4$$','FontName','Times','FontSize',12,...
'interpreter','latex','units','normalized')
text(0.7,0.6,'$$1/\mathrm Q^2$$','FontName','Times','FontSize',12,...
'interpreter','latex','units','normalized')
text(0.7,0.83,'$$1/\mathrm Q^1$$','FontName','Times','FontSize',12,...
'interpreter','latex','units','normalized')
ylabel('$$\left| \mathcal F(\mathrm Q)\right|\,^2$$','interpreter','latex')
text(0.5,-0.1,'$$\mathrm Q R_g$$','interpreter','latex','FontName','Times',...
'FontSize',12,'horizontalalignment','center','units','normalized')
axis square; box on; grid on;
