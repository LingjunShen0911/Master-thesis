function []=plotpoints(PWLModel,Triangles,j)
    if j==0
        i=1;
        figure('visible','off');
    else
        i=j;
    end
    figure(i)
    p = trisurf(Triangles,PWLModel(:,1),PWLModel(:,2),PWLModel(:,3));

    grid on
    set(gca,'LineWidth',2)
    xlim([min(PWLModel(:,1)) max(PWLModel(:,1))])
    xlabel('U_{DS}(V)', 'FontSize', 18, 'FontName', 'Helvetica');
    ylim([min(PWLModel(:,2)) max(PWLModel(:,2))])
    ylabel('U_{GS}(V)', 'FontSize', 18, 'FontName', 'Helvetica');
    zlim([min(PWLModel(:,3)) max(PWLModel(:,3))])
    zlabel('I_{DS}(A)', 'FontSize', 18, 'FontName', 'Helvetica');
    p=patch('Faces',Triangles,'Vertices',PWLModel);

    view(-45,45)
%     view(-35, 20);
     set(p,'FaceColor',[0.95,0.95,0.95],'EdgeColor','black','LineWidth',1)

    if (exist('pics')==0)
        mkdir('pics')
    end
%     set(gca,'FontSize',18, 'FontName','Helvetica')
%     saveas hat neuen Aufruf bei Matlab2015
%     saveas(p,['pics/output.eps'])
%     saveas(p,['pics/output.fig'])
end