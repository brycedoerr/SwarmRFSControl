function movie2gif(mov, gifFile, varargin)
% Movie2Gif ver. 1.0
% ==================
% Matlab movie to GIF Converter.
%
% Syntax: movie2gif(mov, gifFile, prop, value, ...)
% =================================================
% The list of properties is the same like for the command 'imwrite' for the
% file format gif:
%
% 'BackgroundColor' - A scalar integer. This value specifies which index in 
%                     the colormap should be treated as the transparent 
%                     color for the image and is used for certain disposal 
%                     methods in animated GIFs. If X is uint8 or logical, 
%                     then indexing starts at 0. If X is double, then 
%                     indexing starts at 1.
%
% 'Comment' - A string or cell array of strings containing a comment to be
%             added to the image. For a cell array of strings, a carriage 
%             return is added after each row.
%
% 'DelayTime' - A scalar value between 0 and 655 inclusive, that specifies 
%               the delay in seconds before displaying the next image. 
%
% 'DisposalMethod' - One of the following strings, which sets the disposal 
%                    method of an animated GIF: 'leaveInPlace', 
%                    'restoreBG', 'restorePrevious', or 'doNotSpecify'.
%
% 'LoopCount' - A finite integer between 0 and 65535 or the value Inf (the
%               default) which specifies the number of times to repeat the
%               animation. By default, the animation loops continuously. 
%               For a value of 0, the animation will be played once. For a 
%               value of 1, the animation will be played twice, etc. 
%
% 'TransparentColor' - A scalar integer. This value specifies which index 
%                      in the colormap should be treated as the transparent
%                      color for the image. If X is uint8 or logical, then 
%                      indexing starts at 0. If X is double, then indexing
%                      starts at 1
%
% *************************************************************************
% Copyright 2007 Nicolae CINDEA

if (nargin < 2)
  error('Too few input arguments');
end

if (nargin == 2)
    h = waitbar(0, 'Generate GIF file...'); 
    frameNb = size(mov, 2);
    isFirst = true;
    for i = 1:frameNb
        waitbar((i-1)/frameNb, h);
        [RGB, colMap] = frame2im(mov(i)); 
        [IND, map] = aRGB2IND(RGB); 
        if isFirst
            imwrite(IND, map, gifFile, 'gif'); 
            isFirst=false;
        else
            imwrite(IND, map, gifFile, 'gif', 'WriteMode', 'append'); 
        end 
    end
    close(h);
end

if (nargin > 2) 
    h = waitbar(0, 'Generate GIF file...');
    frameNb = size(mov, 2);
    isFirst = true;
    for i = 1:frameNb
        waitbar((i-1)/frameNb, h);
        [RGB, colMap] = frame2im(mov(i)); 
        [IND, map] = aRGB2IND(RGB); 
        if isFirst
            imwrite(IND, map, gifFile, 'gif', varargin{:}); 
            isFirst=false;
        else
            imwrite(IND, map, gifFile, 'gif', 'WriteMode', 'append', ...
                varargin{:}); 
        end 
    end
    close(h);
end