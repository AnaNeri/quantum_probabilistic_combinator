import System.Random (randomRIO)

-- Probabilistic combinator: applies f with probability p, g with probability 1-p
choice :: (a -> a) -> (a -> a) -> Double -> a -> IO a
choice f g p x = do
    prob <- randomRIO (0, 1) :: IO Double
    if prob <= p
        then return (f x)
        else return (g x)

-- Example functions
addOne :: Int -> Int
addOne x = x + 1

subtractOne :: Int -> Int
subtractOne x = x - 1

-- Example usage
main :: IO ()
main = do
    result <- choice addOne subtractOne 0.7 5
    putStrLn $ "Result: " ++ show result