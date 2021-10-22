%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
classdef WebDownload
    properties(Constant)
        HOST='http://cgworkspace.cytogenie.org/';
        PATH='GetDown2/demo';
        KEY_UMAP_EXAMPLES='RunUmapExamples';
        KEY_DEMO='demo';
    end
    
    methods(Static)
        function url=FindUrl(urls)
            app=BasicMap.Global;
            map=app.urlMap;
            
            N=length(urls);
            for i=1:N
                url=urls{i};
                u=java.net.URL(url);
                host=char(u.getHost);
                port=u.getPort;
                if port<=0
                    port=80;
                end
                badKey=['bad:' host ':' num2str(port)];
                if ~map.containsKey(badKey)
                    [ok, issues]=WebDownload.CanReach(host, port);
                    if ok
                        return;
                    end
                    map.put(badKey, badKey);
                end
            end
            if N>0
                url=urls{1};
            else
                url='';
            end
        end
        
        function [ok, issues]=CanReach(host, port, timeout)
            if nargin<3
                timeout=4000;
            end
            issues=java.util.ArrayList;
            ok=edu.stanford.facs.swing.WebDownload.CanReachHost(host, port, timeout, issues);
        end
        
        function url=GetUrlForKey(file, key)
            if nargin<2
                key=WebDownload.KEY_UMAP_EXAMPLES;
                if nargin<1
                    file='';
                end
            end
            app=BasicMap.Global;
            map=app.urlMap;
            url=[];
            if strcmp(key, WebDownload.KEY_UMAP_EXAMPLES)
                if ~map.containsKey(key)
                    url=WebDownload.FindUrl({...
                        [WebDownload.HOST '/run_umap/examples'],...
                        'http://54.39.2.45/run_umap/examples'... %,'https://1drv.ms/u/s!AkbNI8Wap-7_jNJXMO1eLAoecm9ptA?e=Yuvftq/examples'
                        });
                    map.put(key, url);
                else
                    url=map.get(key);
                end 
            elseif strcmp(key, WebDownload.KEY_DEMO)
                if ~map.containsKey(key)
                    url=WebDownload.FindUrl({...
                        [WebDownload.HOST '/GetDown2/demo'],...
                        'http://54.39.2.45/GetDown2/demo'... %,'https://1drv.ms/u/s!AkbNI8Wap-7_jNJXMO1eLAoecm9ptA?e=Yuvftq/examples'
                        });
                    map.put(key, url);
                else
                    url=map.get(key);
                end 
            end
            if ~isempty(url)
                url=[url '/' strrep(file, '\', '/')];
            end
        end
        function url=Url(file, path, host)
            if nargin<3
                host=WebDownload.HOST;
                if nargin<2
                    path=WebDownload.PATH;
                    if nargin<1
                        file='';
                    end
                end
            end
            url=[host path '/' file];
        end
        
        %where only used if Gui.m is present
        function [cancelledByUser, bad, dwl, dlg]=...
                Get(urls, localFiles,  waitWhenDone, allowCancel, where)
            if nargin<5
                where='center';
                if nargin<4
                    allowCancel=true;
                    if nargin<3
                        waitWhenDone=true;
                    end
                end
            end
            try
                dwl=javaObjectEDT('edu.stanford.facs.swing.WebDownload');
            catch ex
                jar=fullfile(fileparts(mfilename('fullpath')), 'webDownload.jar');
                javaaddpath(jar);
                try
                    dwl=javaObjectEDT('edu.stanford.facs.swing.WebDownload');
                catch ex
                    ex.getReport
                    error('Trouble using WebDownload');
                end
            end
            dlg=javaObjectEDT(dwl.dlg);
            try
                Gui.LocateJava(dlg, ...
                    Gui.JWindow(get(0, 'CurrentFigure')), where);
            catch
                WebDownload.LocateJavaOnScreen(dlg, where);
            end
            dwl.waitWhenDone=waitWhenDone;
            dwl.allowCancel=allowCancel;
            cancelledByUser=~dwl.go(urls,localFiles, true);
            bad=dwl.bad;
        end
        
        %where only used if Gui.m is present
        function [ok, cancelledByUser]=GetAndUnzip(url, zipFile, waitWhenDone, ...
                allowCancel, where)
            if nargin<5
                where='center';
                if nargin<4
                    allowCancel=true;
                    if nargin<3
                        waitWhenDone=true;
                    end
                end
            end
            ok=false;
            [cancelledByUser, bad, dwl, dlg]=WebDownload.Get({url}, ...
                {zipFile}, waitWhenDone, allowCancel, where);
            if ~cancelledByUser && bad==0
                dlg.setModal(false);
                dwl.progressBar.setValue(0);
                dlg.setVisible(true);
                dwl.south.setText('Unzipping file now');
                try
                    zipFldr=fileparts(zipFile);
                    if isempty(zipFldr)
                        unzip(zipFile);
                    else
                        unzip(zipFile, zipFldr);
                    end
                    dwl.progressBar.setValue(dwl.progressBar.getMaximum);
                    MatBasics.RunLater(@(h,e)quiet,2);
                    ok=true;
                catch ex
                    ex.getReport
                end
            end
            delete(zipFile);
                    
            function quiet
                dlg.setVisible(false);
            end
        end
        
        function [x,y]=LocateJavaOnScreen(javaComponent, where)
            size=javaComponent.getSize;
            if size.width==0 || size.height==0
                size=javaComponent.getPreferredSize;
            end
            if nargin<2
                where='center';
            end
            p=get(0, 'ScreenSize');
            where=strrep(where,'-', '');
            where=strrep(where,'+','');
            [x, y]=WebDownload.LocateWidthHeight(true, ...
                size.width, size.height, p(1), ...
                p(2), p(3), p(4), where);
            javaComponent.setLocation(x, y);
        end
        
        function [x,y]=LocateWidthHeight(isTop0, width, height, ...
                refX, refY, refWidth, refHeight, where)
            if isempty(where)
                where='center';
            end
            w=String(lower(where));            
            hCenter=~w.contains('east') && ~w.contains('west');
            vCenter=~w.contains('north') && ~w.contains('south');
            centerX=refX+refWidth/2-width/2;
            centerY=refY+refHeight/2-height/2;
            east=w.contains('east');
            top=w.contains('north');
            plusPlus=w.contains('++');
            plus=~plusPlus&&w.contains('+');
            if hCenter && vCenter
                x=centerX;
                y=centerY;
            else
                if width>refWidth
                    if east
                        W=refWidth*.75;
                    else
                        W=refWidth*1.3;
                    end
                else
                    W=width;
                end
                if east
                    if plus
                        h=(refX+refWidth)-W*.4;
                    elseif plusPlus
                        h=(refX+refWidth)-W*.1;
                    else
                        h=(refX+refWidth)-W;
                    end
                else
                    if plus
                        h=refX-(W*.6);
                    elseif plusPlus
                        h=refX-(W*.9);
                    else
                        h=refX;
                    end
                end
                if (top && isTop0) || (~top && ~isTop0)
                    if plus
                        v=refY-height*.25;
                    elseif plusPlus
                        v=refY-height*.9;
                    else
                        v=refY;
                    end
                elseif (top && ~isTop0) || (~top && isTop0)
                    if plus
                        v=(refY+refHeight)-height*.66;
                    elseif plusPlus
                        v=(refY+refHeight)-height*.01;                    
                    else
                        v=(refY+refHeight)-height;
                    end
                end
                if vCenter
                    x=h;
                    y=centerY;
                elseif hCenter
                    x=centerX;
                    y=v;
                else
                    x=h;
                    y=v;
                end
            end
        end
        
        function file=GetZipIfMissing(file, url)
            if ischar(file) && ~exist(file, 'file')
                [fldr,fn]=fileparts(file);
                zipFileName=[fn '.zip'];
                if nargin<2
                    url=fullfile(WebDownload.HOST, WebDownload.PATH);
                end
                url=[url '/' zipFileName];
                zipFile=fullfile(fldr, zipFileName);
                [ok, cancelledByUser]=WebDownload.GetAndUnzip(url, ...
                    zipFile, false, true, 'center');
                if cancelledByUser
                     msg(Html.WrapHr(['Cancelling leaves the required '...
                         'file missing...<br>' BasicMap.Global.smallStart...
                         '<b>' file '</b>' BasicMap.Global.smallEnd]));
                     file=[];
                elseif ~ok
                    msg(Html.WrapHr(['Required file is missing...<br>' ...
                        BasicMap.Global.smallStart...
                        '<b>' file '</b>' BasicMap.Global.smallEnd]));
                    file=[];
                end
            end
        end

    end
end