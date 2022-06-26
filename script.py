import sys
import numpy as np

def gen(n, delta):
  diff = 1
  curr = 0
  xs = [curr]
  store = n
  n -= 1
  while n:
    curr += diff
    xs.append(curr)
    diff = diff + np.random.randint(delta)
    n -= 1
  assert len(xs) == store
  return xs

def check_property(mat):
  n = mat.shape[0]
  max_ = np.zeros((2 * n - 1,))
  min_ = np.full((2 * n - 1,), np.inf)
  for k in range(2 * n - 1):
    for i in range(n):
      j = k - i
      if j < 0 or j >= n:
        continue
      max_[k] = max(max_[k], mat[i, j])
      min_[k] = min(min_[k], mat[i, j])
  for k in range(2 * n - 2):
    print(f'k={k}: max_[{k}]={max_[k]} vs min_[{k+1}]={min_[k+1]}')
    assert max_[k] <= min_[k + 1]

def build_full_matrix(xs):
  n = len(xs)
  mat = np.zeros(shape=(n, n))
  for i in range(n):
    for j in range(n):
      mat[i, j] = xs[i] + xs[j]
  return mat

def main():
  if len(sys.argv) != 3:
    print(f'Usage: python3 {sys.argv[0]} <n:int> <delta:int>')
    sys.exit(-1)
  n = int(sys.argv[1])
  delta = int(sys.argv[2])
  xs = gen(n, delta)
  print(f'xs={xs}')
  mat = build_full_matrix(xs)
  print(f'mat=\n{mat}')
  check_property(mat)

if __name__ == '__main__':
  main()