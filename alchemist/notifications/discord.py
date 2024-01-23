import discord

def completion_alert(token, user_id, message):
    client = discord.Client(intents=discord.Intents.default())

    @client.event
    async def on_ready():
        print(f'Logged in as {client.user} (ID: {client.user.id})')
        user = await client.fetch_user(user_id)
        await user.send(message)
        await client.close()

    client.run(token)
