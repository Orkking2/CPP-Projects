#include <iostream>
#include <deque>
#include <vector>


struct card{
    char id;
    int count = 0; // Default count is 0
};


std::vector<card> IndexCardList(std::deque<card> in, std::vector<char> const &key){ // Sorts to match key
    std::vector<card> CardIndex;
    for(int i = 0; i < sizeof(key); i++){ // For loop here allows for key to change
        CardIndex.push_back(card {key[i]});
    }
    while(!in.empty()){
        card currCard = in.back();
        for(int i = 0; i < sizeof(CardIndex); i++){
            if(currCard.id == CardIndex[i].id){
                CardIndex[i].count++;
                break; // Stops after finding correct id
            }
        }
        in.pop_back();
    }
    return CardIndex;
}

std::vector<char> ReconstructList(std::vector<card> const &IndexedList){
    std::vector<char> out;
    for(card c : IndexedList){
        for(int i = 0; i < c.count; i++){
            out.push_back(c.id);
        }
    }
    return out;
}

void PrintCharList(std::vector<char> in){ // Utility
    for(char c : in){
        std::cout << c << ", ";
    }
    std::cout << "\n";
}

int main()
{
    std::vector<char> key = {'a', '2', '3', '4', '5', '6', '7', '8', '9', '0', 'j', 'q', 'k'}; // Reorder based on frequency increase efficiency

    int l = 200; // Setting random inputs to test sort
    std::deque<card> CardInputList;
    std::vector<char> UnsortedPrint;
    int seed;
    for(int i = 0; i < l; i++){
        CardInputList.emplace_back(card {key[rand() % (key.size() - 1)]});
        UnsortedPrint.emplace_back(CardInputList.back().id);
    }
    PrintCharList(UnsortedPrint);
    //
    PrintCharList(ReconstructList(IndexCardList(CardInputList, key))); // You can unroll these if you want -- just saving lines mostly
}