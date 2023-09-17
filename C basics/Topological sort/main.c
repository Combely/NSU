#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>

typedef struct graphnode {
	int key;
	int indegree;
	int outdegree; 
	struct graphnode** outnodes;
	struct graphnode* next;
}node;

node* graphaddiction_insert(node* head, int prev, int next) {
	node* p = head, *q = head;
	if (!head) {
		//инициализация нового элемента, от которого зависит другой
		head = (node*)malloc(sizeof(node));
		head->outnodes = (node**)malloc(1 * sizeof(node*));
		head->key = prev;
		head->indegree = 0;
		head->outdegree = 1;
		//инициализация нового зависимого элемента
		head->next = (node*)malloc(sizeof(node));
		head->next->outnodes = (node**)malloc(1 * sizeof(node*));
		head->next->next = NULL;
		head->next->key = next;
		head->next->indegree = 1;
		head->next->outdegree = 0;
		//создание абстрактного ребра
		head->outnodes[0] = head->next;
	}
	else {
		while (p->next && p->key != prev) {
			p = p->next;
		}
		while (q->next && q->key != next) {
			q = q->next;
		}
		if (q->key == next) {
			if (p->key == prev) { //случай, когда элементы и исходящей, и принимающей вершин ребра есть в списке
				//проверка, была ли уже запись с такой парой чисел
				int flag = 0;
				for (int i = 0; i < p->outdegree; i++)
					if (q->key == p->outnodes[i]->key)
						flag = 1;
				//создать новую зависимость, если ее еще не было
				if (!flag) {
					p->outdegree++;
					p->outnodes = (node**)realloc(p->outnodes, p->outdegree * sizeof(node*));
					p->outnodes[p->outdegree - 1] = q;
					q->indegree++;
				}
			}
			else if (!p->next) { //случай, когда элемент принимающей вершины ребра есть в списке, а исходящей - нет
				p->next = (node*)malloc(sizeof(node));
				p->next->outnodes = (node**)malloc(1 * sizeof(node*));
				p->next->key = prev;
				p->next->indegree = 0;
				p->next->outdegree = 1;
				p->next->outnodes[p->next->outdegree - 1] = q;
				p->next->next = NULL;
				q->indegree++;
			}
		}
		else if (!q->next) {
			if (p->key == prev) {  //случай, когда элемент исходящей вершины ребра есть в списке, а принимающей - нет
				//инициализация нового элемента, который зависит от существующего
				q->next = (node*)malloc(sizeof(node));
				q->next->outnodes = (node**)malloc(1 * sizeof(node*));
				q->next->key = next;
				q->next->indegree = 1;
				q->next->outdegree = 0;
				q->next->next = NULL;
				//обновление существующего элемента
				p->outdegree++;
				p->outnodes = (node**)realloc(p->outnodes, p->outdegree * sizeof(node*));
				p->outnodes[p->outdegree - 1] = q->next;
			}
			else if (!p->next) {  //случай, когда элементов и исходящей, и принимающей вершин ребра нет в списке
				p->next = (node*)malloc(sizeof(node));
				p->next->key = prev;
				p->next->indegree = 0;
				p->next->outdegree = 1;
				p->next->outnodes = (node**)malloc(1 * sizeof(node*));
				q = p->next;
				q->next = (node*)malloc(sizeof(node));
				q->next->key = next;
				q->next->indegree = 1;
				q->next->outdegree = 0;
				q->next->outnodes = (node**)malloc(1 * sizeof(node*));
				q->next->next = NULL;
				p->next->outnodes[p->next->outdegree - 1] = q->next;
			}
		}
	}	
	return head;
}

node* singlegnode_insert(node* head, int num) {
	node* p = head;
	while (p->next) {
		p = p->next;
	}
	p->next = (node*)malloc(sizeof(node));
	p->next->key = num;
	p->next->indegree = 0;
	p->next->outdegree = 0;
	p->next->outnodes = NULL;
	p->next->next = NULL;
	return head;
}

node* kahntpsort(node* head) { 
	node* newhead = NULL, * curr = head, * minfreenode = NULL, * prev = NULL, * prevformin = NULL;
	while (curr->next) {
		if (!curr->indegree) {
			if (minfreenode) {
				if (curr->key < minfreenode->key) {
					prevformin = prev;
					minfreenode = curr;
				}
			}
			else {
				prevformin = prev;
				minfreenode = curr;
			}
		}
		prev = curr;
		curr = curr->next;
	}
	if (!curr->indegree) {
		if (minfreenode) {
			if (curr->key < minfreenode->key) {
				prevformin = prev;
				minfreenode = curr;
			}
		}
		else {
			prevformin = prev;
			minfreenode = curr;
		}
	}
	if (minfreenode) {
		//удаление всех выходящих ребер и уменьшение числа входивших в вершины из удаляемой ребер
		if (minfreenode->outnodes) { //не будет выполняться, если minfreenode - это вершина без зависимостей
			for (int i = 0; i < minfreenode->outdegree; i++) {
				minfreenode->outnodes[i]->indegree--;
			}
			free(minfreenode->outnodes);
			minfreenode->outnodes = NULL;
		}
		//удаление элемента из предыдущего списка
		if (minfreenode == head) {
			head = minfreenode->next;
		}
		else if (!minfreenode->next) {
			prevformin->next = NULL;
		}
		else {
			prevformin->next = minfreenode->next;
		}
		newhead = minfreenode;
		if (head) {
			//условие на head нужно, чтобы не перепутать завершающий случай отсутствия элементов в изначальном
			//списке со случаем отсутствия элементов с indegree == 0
			newhead->next = kahntpsort(head);
			if (!newhead->next)
				newhead = NULL;
		}
		else {
			newhead->next = NULL;
		}
	}
	//Старый список, на который указывает head, должен с каждым вызовом функции kahntpsort уменьшаться на 1 элемент.
	//Возвращать же функция должна минимальный элемент, в который не входят никакие вершины.
	//Новый список должен выстраиваться по первому вызову newhead и следующим за ним элементам.
	//Если список head не нулевой, но в нем не нашлось вершины с 0 входящих ребер, то рекурсию не вызывать.
	//Эта идея основывается на том, что в ациклическом графе при удалении каждого элемента списка вплоть до
	//head == NULL будет хоть одна вершина с 0 входящих ребер.
	//Как из вызовов функции понять, есть ли циклы. В ациклическом графе до тех пор, пока изначальный список не станет
	//равным нулю, будет возвращаться минимальный элемент из тех, которые имеют indegree == 0. Если же имеется
	//внутренний цикл в графе, то весь построенный новый список (если он есть) рекурсивно занулится.
	return newhead;
}

void clearoff(node* head) {
	node* p;
	while (head->next) {
		p = head;
		head = head->next;
		free(p);
	}
	free(head);
}

int main() {
	FILE* in = fopen("input.txt", "r"), * out = fopen("output.txt", "w");
	int N, M, firsttheme, sectheme, *addictionpart;
	node* head = NULL, *sortedgraphhead;
	fscanf(in, "%d%d", &N, &M);
	addictionpart = (int*)malloc(N * sizeof(int));
	for (int i = 0; i < N; i++) {
		addictionpart[i] = 0;
	}
	for (int i = 0; i < M; i++) {
		fscanf(in, "%d%d", &firsttheme, &sectheme);
		addictionpart[firsttheme - 1] = 1;
		addictionpart[sectheme - 1] = 1;
		head = graphaddiction_insert(head, firsttheme, sectheme);
	}
	for (int i = 0; i < N; i++) {
		if (addictionpart[i] == 0)
			head = singlegnode_insert(head, i + 1);
	}
	sortedgraphhead = kahntpsort(head);
	/*node* p = sortedgraphhead;
	while (p->next) {
		printf("%d ", p->key);
		p = p->next;
	}
	printf("%d ", p->key);*/
	if (sortedgraphhead == NULL) { //если kahntpsort вернул NULL, значит в head есть циклы
		fprintf(out, "bad course");
	}
	else {
		node* p = sortedgraphhead;
		while (p->next) {
			fprintf(out, "%d ", p->key);
			p = p->next;
		}
		fprintf(out, "%d ", p->key);
	}
	if (sortedgraphhead)
		clearoff(sortedgraphhead);
	free(addictionpart);
	fclose(in);
	fclose(out);
	return 0;
}