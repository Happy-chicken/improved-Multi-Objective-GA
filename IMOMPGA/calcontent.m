function content = calcontent(time,fuzzy_mat)
%time完成时间
%fuzzy_mat 模糊隶属度矩阵

%content满意度
[row,~]=size(fuzzy_mat);
content = zeros(1,6);
for i = 1:row
    if (time(i)<=fuzzy_mat(i,1)) || (time(i)>fuzzy_mat(i,4))
        content(i) = 0;
    elseif (time(i)>fuzzy_mat(i,1)) && (time(i)<=fuzzy_mat(i,2))
        content(i) = (time(i) - fuzzy_mat(i,1))/(fuzzy_mat(i,2)-fuzzy_mat(i,1));
    elseif (time(i)>fuzzy_mat(i,2)) && (time(i)<=fuzzy_mat(i,3))
        content(i) = 1;
    elseif (time(i)>fuzzy_mat(i,3)) && (time(i)<=fuzzy_mat(i,4))
        content(i) = (fuzzy_mat(i,4)-time(i))/(fuzzy_mat(i,4)-fuzzy_mat(i,3));
    end
end
end

