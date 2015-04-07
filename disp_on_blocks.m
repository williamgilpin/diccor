function [bx,by,corsh] = disp_on_blocks(im1,im2,blocksize,subpix)

% DISP_ON_BLOCKS calculates x- and y-displacements for each distinct block 
%	taken from the image im1. The size of the blocks is determind by the 
%	parameter blocksize. For each block in the image im1 the corresponding
%	block in the image im2 is found (i.e., the block at the same position),
%	and the cross-correlation function between the two blocks is formed.
%	Coordinates of the peak of the cross-correlation function constitute
%	the displacement vector of the block. (Displacements go from im1 to im2.) 
%	The parameter subpix determines whether the resulting values of displacements 
%	are integers or not; if subpix > 0, the values of displacements are 
%	non-integers, otherwise they are integers.

%	bx is a matrix of x-displacements, by a matrix of y-displacements, and
%	corsh the fft-shifted cross-correlation by block between im1 and im2.

%   Xavier Trepat 07-07.
%	Iva Marija Tolic-Norrelykke 03-21-01

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

if nargin < 4,
   subpix	= 	0;
end;
% SUBTRACT THE MEAN OF EACH BLOCK
im1 = blkproc(im1,[blocksize blocksize],'x - mean2(x)');
im2 = blkproc(im2,[blocksize blocksize],'x - mean2(x)');
% myfun = @(x) x.data-mean2(x.data);
% im1 = blockproc(im1,[blocksize blocksize],myfun);
% im2 = blockproc(im2,[blocksize blocksize],myfun);
% % % MULTIPLY BY A 2D HANNING WINDOW
% im1 = blkproc(im1,[blocksize blocksize],'xhanning2D(x)');
% im2 = blkproc(im2,[blocksize blocksize],'xhanning2D(x)');

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% CALCULATE THE CROSS-CORRELATION FUNCTION

autocor_whole1 = blkproc(im1,[blocksize blocksize],'ifft2((fft2(x)).*conj(fft2(x)))');
autocor_whole2 = blkproc(im2,[blocksize blocksize],'ifft2((fft2(x)).*conj(fft2(x)))');
autocor_max1   = blkproc(autocor_whole1,[blocksize blocksize],'repmat(x(1,1),size(x))');
autocor_max2   = blkproc(autocor_whole2,[blocksize blocksize],'repmat(x(1,1),size(x))');
cor_whole1     = blkproc(im1,[blocksize blocksize],'conj(fft2(x))');
cor_whole2     = blkproc(im2,[blocksize blocksize],'fft2(x)');
cor_whole      = cor_whole2 .* cor_whole1;  % takes the norm of the fourier transform at each point
cor_whole		= blkproc(cor_whole,[blocksize blocksize],'ifft2(x)');
cor_whole      = real(cor_whole) ./ sqrt(autocor_max1.*autocor_max2);
cor_whole(find(~isfinite(cor_whole))) = zeros(size(find(~isfinite(cor_whole))));
corsh          = blkproc(cor_whole,[blocksize blocksize],'fftshift(x)');
%imagesc(corsh);
%whos corsh
%plot(corsh(100,:))
%user_entry = input('prompt') 
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% FIND THE PEAK OF THE CROSS-CORRELATION FUNCTION

 %if subpix > 0,
%  bx = blkproc(corsh,[blocksize blocksize],'x_cntr(x,2)');
%  by = blkproc(corsh,[blocksize blocksize],'y_cntr(x,2)');
	bx=blkproc(corsh,[blocksize blocksize],'center_x_1d(x)');
	by=blkproc(corsh,[blocksize blocksize],'center_y_1d(x)');
%whos cliffxout
%user_entry = input('prompt') 
    %else   
%	bx				= blkproc(corsh,[blocksize blocksize],'max_cor_bx(x)');
%	by				= blkproc(corsh,[blocksize blocksize],'max_cor_by(x)');
%end; %(if subpix > 0)