function swap(list, x, y) {
    var b = list[y];
    list[y] = list[x];
    list[x] = b;
}

function partition(list, low, high){
    var pivot = list[high];
    var i = low - 1;

    for (let j = low; j <= high ; j++) {
        if(list[j] <= pivot){
            i++;
            swap(list, i, j);
        }
    }
    swap(list, i+1, high);
    return i+1;
}