# 재료역학 기말 프로젝트 파일
재료역학 과제입니다.
# 실행 방법
`cmake`를 이용해서 C++ 프로젝트를 빌드합니다. 입력은 input/input.txt에 넣어두십시오. 준비가 완료되었으면 다음과 같은 순서로 실행합니다.
```
cd ./build
cmake ../
make
./main
```
이 때 프로젝트의 루트 디렉토리에 `data.txt`가 있으면 빌드에 성공한 것입니다.
# 그래프 그리기
`requirements.txt`에 있는 파이썬 라이브러리를 설치합니다.
```
python drawSFD.py
python drawBMD.py
```
를 통해서 그래프를 그립니다. 만들어진 그래프는 `./output/`에 pdf 형태로 저장되어 있습니다.

