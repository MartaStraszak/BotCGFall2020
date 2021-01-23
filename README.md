# BotCGFall2020
My bot for the [CodinGame Fall Challenge 2020](https://www.codingame.com/contests/fall-challenge-2020). Reached the Legend rank and placed 31th out of 7000 participants. See the [final ranking](https://www.codingame.com/contests/fall-challenge-2020/leaderboard/global) and [my profile](https://www.codingame.com/profile/925c8519d0f84c9dba6c86560f439f068401402).

# Basic Insights
## First learn then cast/brew
One can observe very early on that a simple, yet quite effective strategy is to start by buying several (something around `~10` seems optimal) spells and then start casting spells and brewing potions without buying any new spells anymore. While I have tested more sophisticated strategies and actually have written a lot of code to implement seemingly smarter ideas, I can confidently say that it was enough to stick with the basic thing to reach Legend.

## Ignore the Opponent

Another conceptual simplification which seems counterintuitive at first but brings a lot of clarity and allows to focus on simpler goals is to treat each move as a separate instance of a single-player game. In fact, with the exception of a few absolutely top-scoring bots, this was the usual assumption among high-ranked player. The rationale behind this assumption is that taking into account the other player and the impact of his actions, is simply too complex and actually not worth it. We will use different techniques to ensure that our bot is not too greedy and actually tries to brew potions when possible, without delaying rewards too much.

In practice, "ignoring the opponent" simply means to find a, say, 8 (or 10 or even longer) consecutive moves from the current position to maximize some kind of reward gained during these 8 moves. During this 8 move we treat the opponent as "frozen" meaning that he will not cast any spells and importantly not brew any potions (to potentially steal rewards from us).

# The Actual Algorithm

## Finding a Move via Optimization
From the above basic insights, the problem we need to solve is given a state `S` in the game consisting of:

 - The counts of each ingredient `[c1, c2, c3, c4]`,
 - The potions available for brew along with the rewards for brewing them,
 - The set of spells we have learned,
we would like to generate a sequence of `K` moves (in my final code `K=14` but even `K=10` or `K=8` are still fine) `[m1, m2, m3, ..., mK]` that maximizes the total reward of such a "path". The move we actually make at this state is then  of course, simply `m1`. The moves `m_i` are either `cast`, `brew` or `rest` and the reward is the sum of the following:
 - The sum of points earned for brewing potions during these `K` moves
 - The "inventory bonus" which is the sum `w1*c1'+w2*c2'+w3*c3'+w4*c4'` where `[c1', c2', c3', c4']`is the final state of the inventory (counts of ingredient after the `K` moves).

By maximizing the score of such a path the bot intuitively picks the best possible strategy not just "for now" but for the "nearest future" (in this case `K` moves), which of course can be much more sophisticated then greedily going for the quickest brew or so.

One issue that such a bot would suffer from is that it is actually too relaxed about "ignoring the opponent". This means that if there are two paths that give roughly equal inventory bonus, both brewing a single (the same) potion, but one of them brews during the `3rd` move and the second one brews during the `11th` move, then they will have pretty much the same scores. However the first path is clearly better as the opponent is not so likely to steal the potion from us then. For this reason, as commonly done in Reinforcement Learning and AI in general, we introduce a "discount factor" `0<gamma<1` (think `gamma=0.97`) and "weigh" the reward gained at step `k` with `gamma^k`, thus in other words, the later the reward is gained the less does it contribute to the score of a path.

## Finding the Best Path -- Beam Search

Having the above formulation saying that we need to find a path that maximizes its score, one can try to apply some standard optimization techniques to solve this problem. What seemed to work best for me in this case was beam search. You essentially search the best path using bruteforce (try each possible option at each intermediate state) but, after each new step, prune the set of "partial paths" so that some notion of "partial score" is best possible. This is of course a heuristic method and thus also some engineering comes into how to design this "partial score" function, so feel free to take a closer look at the source code in C++.

To make the code run in the allowed time for reasonable `K`'s one has to do a lot of optimizations, such as precomputing transitions between inventory states when casting a particular spell or brewing a potion and using bitsets to represent arrays of bools whenever possible (and use bit operations on them to speed things up).


# More Sophisticated Ideas

Here is a list of ideas that I have tried but they did not make the bot much stronger (to the extent I could notice that in local and on-server tests):

 1. One can assign some scores to spells that are available for buy in order to estimate their overall usefullness. Since the set of spells is static one can use some historical data (for instance from the self-play of the bot against himself) to precompute which spells are good and which are not. Empirically, either my method for assigning scores, or the whole idea did not work that well...
 2. Adding the "learn spell" action as one of the available moves in the beam search is something you can do, and as long as you allow this move only in the several initial steps (`2` is fine) on the path, this does not seem to make the search much slower compared to the standard one. However, the gain from doing that did not seem so great. In fact, even though it ended up being part of my finall submission I don't think it was a game-changer. In particular, I think the bot would reach Legend even without this improvement.
 3. Actually looking at what the opponent does and taking this into account during the beam search was the next thing I wanted to implement, but due to lack of time I didn't manage to do that. This seems to be something that could help a lot, but also would require significant work to do properly.

