function mypdf(fname,r,s)

if nargin < 2,
    r = .71; % height/width ratio
end;
if nargin < 3, 
    s = 1; % scaling of font size
end;


set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperSize',s*[13.2,r*13.2]);
set(gcf,'PaperPosition',s*[.1,.1,13,r*13]);
print(gcf,'-dpdf', fname);
%print(gcf,'-djpeg',fname);
print(gcf,'-depsc2',fname);



% if strfind(fname,'Loss')
%     set(legend(gca),'FontSize',12)
% end
% 
% 
% opts = struct('bounds','loose','color','cmyk','FontSize',1.2*s);
% exportfig(gcf,fname,opts);
% % print(gcf,'-depsc',fname);
