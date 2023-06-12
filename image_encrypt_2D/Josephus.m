%% 约瑟夫（josephus）问题
function result = Josephus(num,start,step)   %num 总数  start开始点  step 步长
    node = 1:num;
    now = start;  %now记录当前在哪个人
    call = 1;     %call记录当前叫号
    result = [];  %存放结果序列
    while length(node) ~= 0
        if mod(call,step)==0
            result = [result,now];
            call = 1;
            now_1 = now;      %暂存当前位置 进行删除
            if now==max(node)
                now = node(1);
            else
                now = node(find(node==now)+1);
            end
            node(find(node==now_1)) = [];
        else
            call = call+1;
            if now==max(node)
                now = node(1);
            else
                now = node(find(node==now)+1);
            end
        end
    end
end