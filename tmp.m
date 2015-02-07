 I = imread('circuit.tif');
       BW1 = edge(I,'prewitt');
       BW2 = edge(I,'canny');
       figure, imshow(BW1)
       figure, imshow(BW2)