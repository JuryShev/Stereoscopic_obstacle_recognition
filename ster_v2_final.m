%% Загрузка стереопары
Left=(imread('left_1.60.jpg'));
Righ=(imread('righ_1.60.jpg'));
% избавление от искжений (от дисторсии)
Left=undistortImage(Left, stereoParams2.CameraParameters1);
Righ=undistortImage(Righ, stereoParams2.CameraParameters2);
Left_gray=rgb2gray(Left);
Righ_gray=rgb2gray(Righ);
% выравнивание изображений
[L,R] = rectifyStereoImages(Left_gray,Righ_gray,stereoParams2);

%% 
%Интерполяция изображения 
L2=imresize(L,0.5,'nearest');
R2=imresize(R,0.5,'nearest');
% % %Медианная фильтрация
L2=medfilt2(L2,[7, 7]);
R2=medfilt2(R2, [7, 7]);
%% Расчет карты глубины изображения
disparityRange = [0 80];
disparityMap = disparity(L2,R2,'BlockSize',...
   9,'DisparityRange',disparityRange );
%% 
disparityMap=uint8(disparityMap);
disparityMap_G=mat2gray(disparityMap,disparityRange);
disparityMap_G2 = imfill(disparityMap_G,'holes');
disparityMap_G2=medfilt2(disparityMap_G2,[5, 5]);
figure, imshow(disparityMap_G2)
%% Нахождение углов в интересующей области
points1 = detectHarrisFeatures(L, 'ROI', [991, 367, 1239, 1003 ],'MinQuality', 0.3);
points2 = detectHarrisFeatures(R, 'ROI',[991, 367, 1239, 1003 ],'MinQuality', 0.3);
%% Извлечение функции окрестности.
[features1,valid_points1] = extractFeatures(L,points1);
[features2,valid_points2] = extractFeatures(R,points2);
%% сопастовление ключевых точек
indexPairs = matchFeatures(features1,features2);
%%
matchedPoints1 = valid_points1(indexPairs(:,1),:);
matchedPoints2 = valid_points2(indexPairs(:,2),:);
%% Вывод координат точек на изображении
k1=matchedPoints1.Location;
k2=matchedPoints2.Location;
%% триангуляция
point3d = triangulate(k1, k2, stereoParams2);
k1(:,3:4)=20;
k2(:,3:4)=20;
%% сортировка sort_k1

nn=1;
for n=1:length(k1)
    if (k1(n,1)>=991&&k1(n,1)<=2.23e+03)&&(k1(n,2)>=367&&k1(n,2)<1.37e+03)&&point3d(n,3)>0
        k1_sort(nn,1)=k1(nn,1);
        k1_sort(nn,2)=k1(nn,2);
        k1_sort(nn,3)=n;
        k1_sort(nn,4)= point3d(n,3);
        nn=nn+1;
    end  
end

%% Сегментация объекта
k1_sort=sortrows(k1_sort,4);
I=disparityMap_G2;
I=imresize(I,2,'nearest');

%% Выделение интересующей области по яркости пекселя
seedpointR = round(k1_sort(1,2));
seedpointC = round(k1_sort(1,1));
%% Сегментация
W = graydiffweight(I, seedpointC, seedpointR,'GrayDifferenceCutoff',40);
figure, imshow(log(W),[])
%% Перевод в бинарное изображение
thresh = 0.001;
BW = imsegfmm(W, seedpointC, seedpointR, thresh);
figure, imshow(BW)
title('Segmented Image')
%% Морфологическое преобразование
se = strel('disk',4);
erodedI = imopen(BW,se);

%% Вычисление границ области сегментации

Box=regionprops(erodedI,'BoundingBox');
for k = 1 : length(Box)
  thisBB = Box(k).BoundingBox;
  M=rectangle('Position', [thisBB(1),thisBB(2),thisBB(3),thisBB(4)],...
  'EdgeColor','r','LineWidth',2 );

end
%% вычисление дистанции и  масштаба пикселей
distanseInMeters=(point3d(k1_sort(1,3),3)/1000)*1.2;
scale=abs(point3d(k1_sort(1,3),1)-point3d(k1_sort(2,3),1))/abs(k1(k1_sort(1,3), 1)-k1(k1_sort(2,3), 1));% масштаб в 1 пиксель 1.63
%% Выод результата
for k = 1 : length(Box)
thisBB = Box(k).BoundingBox;
W=thisBB(1,3)*scale/1000;%длинна области
H=thisBB(1,4)*scale/1000;%высота области
    if W*H>0.01
    distanceAsString = sprintf('Distanse = %0.2fmeters W = %0.2fmeters H = %0.2fmeters', distanseInMeters, W,H);
    L = insertObjectAnnotation(L,'rectangle',thisBB ,distanceAsString,'FontSize', 40, 'LineWidth', 10);
    end
end
figure, imshow(L)
